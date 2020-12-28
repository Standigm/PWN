We also use a protein-protein interaction network from the [BioGRID](https://thebiogrid.org/) for additional
experiments.
To reproduce `adj.npz` and `ensp.txt`, first install some Python requirements:

```bash
> pip install "pwn[data] @ git+https://github.com/Standigm/PWN"
```

Then run the following Python script:

```python
import io
from collections import defaultdict
from zipfile import ZipFile

import requests
import scipy.sparse as sp
import pandas as pd
from tqdm import tqdm


idmap = (
    pd.read_csv(
        'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.203/BIOGRID-IDENTIFIERS-4.4.203.tab.zip',
        skiprows=28,
        sep='\t',
    )
    .query('IDENTIFIER_TYPE == "ENSEMBL" and ORGANISM_OFFICIAL_NAME == "Homo sapiens"')
    .rename(columns={'BIOGRID_ID': 'biogrid', 'IDENTIFIER_VALUE': 'ensg'})[['biogrid', 'ensg']]
    .set_index('biogrid')
)


df = pd.read_csv(
    ZipFile(
        io.BytesIO(
            requests.get(
                'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.203/BIOGRID-ORGANISM-4.4.203.tab3.zip'
            ).content
        )
    ).open('BIOGRID-ORGANISM-Homo_sapiens-4.4.203.tab3.txt'),
    sep='\t',
    usecols=['BioGRID ID Interactor A', 'BioGRID ID Interactor B', 'Experimental System Type', 'Score'],
)
edges = (
    df.assign(Score=df.Score.map(lambda x: float(x if x != '-' else 'nan')))
    .query('`Experimental System Type` == "physical"')
    .set_index('BioGRID ID Interactor A')
    .join(idmap, how='inner')
    .rename(columns={'ensg': 'src'})
    .set_index('BioGRID ID Interactor B')
    .join(idmap, how='inner')
    .rename(columns={'ensg': 'dst'})[['src', 'dst']]
)
nodes = sorted(set(edges['src']) | set(edges['dst']))

node2idx = {node: i for i, node in enumerate(nodes)}
adj = sp.dok_matrix((len(nodes), len(nodes)), dtype=bool)
for _, src, dst in tqdm(edges.itertuples(), total=len(edges)):
    i, j = node2idx[src], node2idx[dst]
    adj[i, j] = adj[j, i] = True

labels = sp.csgraph.connected_components(adj)[1]
clustersize = defaultdict(int)
for l in labels:
    clustersize[l] += 1
cluster = max(clustersize, key=lambda x: clustersize[x])
select = labels == cluster

pd.DataFrame([node for cond, node in zip(select, nodes) if cond]).to_csv('ensg.txt', index=False, header=None)
sp.save_npz('adj.npz', adj.tocsr()[select].tocsc()[:, select].tocoo())
```

## License

The data is originally provided by the Mike Tyers Lab and redistributed under a
[MIT License](https://biogrid-downloads.nyc3.digitaloceanspaces.com/LICENSE.txt).
