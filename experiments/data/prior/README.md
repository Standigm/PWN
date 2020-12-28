We use a prior knowledge from the [Cancer Gene Census](https://cancer.sanger.ac.uk/census) (GRCh38, v92).
To reproduce `ensg.txt`, first download `cancer_gene_census.csv` from [COSMIC](https://cancer.sanger.ac.uk/cosmic).
Then, install some Python requirements:

```bash
> pip install "pwn[data] @ git+https://github.com/Standigm/PWN"
```

Finally, run the following Python script:

```python
import pandas as pd
from tqdm import tqdm

prior = []
df = pd.read_csv('cancer_gene_census.csv', usecols=['Synonyms'], dtype=str).dropna()
for _, synonyms in tqdm(df.itertuples(), total=len(df)):
    ensg = [gene.split('.')[0] for gene in synonyms.split(',') if gene.startswith('ENSG')]
    assert len(ensg) == 1
    prior.append(ensg[0])

pd.DataFrame({'ENSG': prior}).to_csv('ensg.txt', index=False, header=None)
```

## License

The data is originally provided by the Wellcome Sanger Institute.
See their [legal terms and conditions](https://www.sanger.ac.uk/legal/) for more details.
