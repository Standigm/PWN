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


def generate_query(idtype: str, identifiers: Iterable[str]) -> str:
    if idtype == 'ensp':
        idtype = 'ensembl_peptide_id'
    elif idtype == 'uniprot':
        idtype = 'uniprot_gn_id'
    else:
        raise NotImplementedError
    API = 'http://www.ensembl.org/biomart/martservice?query='
    QUERY = f"""<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
    <Query virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6">
        <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
            <Filter name = "{idtype}" value = "{{}}"/>
            <Attribute name = "ensembl_gene_id" />
            <Attribute name = "{idtype}" />
        </Dataset>
    </Query>""".replace(
        '\n', ''
    )
    return re.sub(' +', ' ', f'{API}{QUERY}').replace(' ', '%20').format(','.join(identifiers))


for ppi, idtype in {'string': 'ensp', 'iid': 'uniprot'}.items():
    identifiers = set(pd.read_csv(f'../ppi/{ppi}/{idtype}.txt', header=None)[0])
    df = pd.concat(
        [pd.read_csv(generate_query(idtype, ids), sep='\t', header=None) for ids in tqdm(list(grouper(identifiers, 256)))]
    ).drop_duplicates()
    df.to_csv(f'{ppi}.tsv', index=False, header=None, sep='\t')
```

## License

The data is originally provided by the European Bioinformatics Institute at the European Molecular Biology Laboratory.
See the [disclaimer](https://www.ensembl.org/info/about/legal/disclaimer.html) for more details.
