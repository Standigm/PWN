Since the gene ID types are different, we need a ID converter.
To reproduce `idmap.tsv`, first install some Python requirements:

```bash
> pip install "pwn[data] @ git+https://github.com/Standigm/PWN"
```

Then run the following Python script to fetch data from the [Ensembl](https://www.ensembl.org/index.html):

```python
import re
from collections.abc import Iterable

import pandas as pd
from tqdm import tqdm


def grouper(iterable: Iterable, n: int) -> Iterable[list]:
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3) --> ABC DEF G
    from itertools import zip_longest

    args = [iter(iterable)] * n
    return map(lambda xs: list(filter(lambda x: x is not None, xs)), zip_longest(*args))


def generate_query(ensp: Iterable[str]) -> str:
    API = 'http://www.ensembl.org/biomart/martservice?query='
    QUERY = """<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
    <Query virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6">
        <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
            <Filter name = "ensembl_peptide_id" value = "{}"/>
            <Attribute name = "hgnc_symbol" />
            <Attribute name = "ensembl_gene_id" />
            <Attribute name = "ensembl_peptide_id" />
        </Dataset>
    </Query>""".replace(
        '\n', ''
    )
    return re.sub(' +', ' ', f'{API}{QUERY}').replace(' ', '%20').format(','.join(ensp))


ensp = pd.concat([pd.read_csv(f'../ppi/{source}/ensp.txt', header=None) for source in ['string']], ignore_index=True)
df = pd.concat(
    [pd.read_csv(generate_query(ensp), sep='\t', header=None) for ensp in tqdm(list(grouper(set(ensp[0]), 256)))]
).drop_duplicates()
df.to_csv('idmap.tsv', index=False, header=None, sep='\t')
```

## License

The data is originally provided by the European Bioinformatics Institute at the European Molecular Biology Laboratory.
See the [disclaimer](https://www.ensembl.org/info/about/legal/disclaimer.html) for more details.
