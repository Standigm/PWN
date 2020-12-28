We use a omics data from [The Cancer Genome Atlas](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga).
To reproduce `pvalue.csv` and `log2fc.csv`, first download `gdc-client` from the [Genomic Data Commons](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
and execute following commands to download raw data:

```bash
> gdc-client download -m metadata/gdc_manifest_20201013_122018.txt -d rawdata  # takes almost 12 hours
```

Then, install some Python requirements:

```bash
> pip install "pwn[data] @ git+https://github.com/Standigm/PWN"
```

Finally, run the following Python script (which takes almost 30 minutes):

```python
import numpy as np
import scipy.stats as ss
import pandas as pd
from tqdm import tqdm

tqdm.pandas()


metadata = (
    pd.read_csv('metadata/gdc_sample_sheet.2020-10-13.tsv', sep='\t')
    .query('`Sample Type` == "Primary Tumor" or `Sample Type` == "Solid Tissue Normal"')
)
metadata = pd.DataFrame.from_dict(
    {
        'id': metadata['Case ID'],
        'filepath': metadata['File ID'] + '/' + metadata['File Name'],
        'disease': metadata['Project ID'].astype('category'),
        'control': metadata['Sample Type'] == 'Solid Tissue Normal',
    }
)
metadata = pd.merge(metadata.id[metadata.control], metadata.id[~metadata.control]).drop_duplicates().merge(metadata)
metadata = metadata.groupby('id').count().query('disease == 2').reset_index()[['id']].merge(metadata)
metadata = metadata.groupby('disease').count().query('id >= 30').reset_index()[['disease']].merge(metadata)

dfs = []
for row in tqdm(metadata.itertuples(), total=len(metadata)):
    df = pd.read_csv('rawdata/' + row.filepath, sep='\t', names=['ensg', 'value'])
    df['ensg'] = list(map(lambda x: x.split('.')[0], df['ensg']))
    df['patient'] = row.id
    df['disease'] = row.disease
    df['control'] = row.control
    df['value'] = np.log2(df['value'] + 1.0)
    dfs.append(df)
df = pd.concat(dfs, ignore_index=True)

df = df.pivot(index=['ensg', 'patient', 'disease'], columns='control', values='value')
df = df.assign(log2fc=df[False] - df[True])['log2fc'].reset_index().groupby(['disease', 'ensg'])
(
    df.progress_apply(lambda df: ss.ttest_1samp(df.log2fc, 0.0).pvalue)
    .to_frame('pvalue')
    .reset_index()
    .pivot(index='ensg', columns='disease', values='pvalue')
    .to_csv('pvalue.csv')
)
(
    df.agg({'log2fc': 'mean'})['log2fc']
    .to_frame('log2fc')
    .reset_index()
    .pivot(index='ensg', columns='disease', values='log2fc')
    .to_csv('log2fc.csv')
)
```

## License

The data is originally provided by the Center for Cancer Genomics at the National Cancer Institute.
See the [data access policies](https://gdc.cancer.gov/access-data/data-access-policies) for more details.
