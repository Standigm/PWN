We also use a protein-protein interaction network from the [IID](http://iid.ophid.utoronto.ca) for additional
experiments.
To reproduce `adj.npz` and `uniprot.txt`, first install some Python requirements:

```bash
> pip install "pwn[data] @ git+https://github.com/Standigm/PWN"
```

Then run the following Python script:

```python
from collections import defaultdict

import scipy.sparse as sp
import pandas as pd
from tqdm import tqdm


df = (
    pd.read_csv(
        'http://iid.ophid.utoronto.ca/static/download/human_annotated_PPIs.txt.gz',
        sep='\t',
        usecols=['uniprot1', 'uniprot2', 'evidence_type'],
    )
)
edges = (
    df
    .query('evidence_type.str.startswith("exp")')
    .rename(columns={'uniprot1': 'src', 'uniprot2': 'dst'})[['src', 'dst']]
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

pd.DataFrame([node for cond, node in zip(select, nodes) if cond]).to_csv('uniprot.txt', index=False, header=None)
sp.save_npz('adj.npz', adj.tocsr()[select].tocsc()[:, select].tocoo())
```
