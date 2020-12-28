We use a protein-protein interaction network from the [STRING](https://string-db.org).
To reproduce `adj.npz` and `ensp.txt`, first install some Python requirements:

```bash
> pip install "pwn[data] @ git+https://github.com/Standigm/PWN"
```

Then run the following Python script:

```python
from collections import defaultdict

import scipy.sparse as sp
import pandas as pd
from tqdm import tqdm


edges = (
    pd.read_csv(
        'https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz',
        sep=' ',
        usecols=['protein1', 'protein2', 'combined_score', 'experiments', 'database'],
    )
    .query('combined_score > 700 & ((experiments > 0) | (database > 0))')
    .rename(columns={'protein1': 'src', 'protein2': 'dst'})[['src', 'dst']]
    .apply(lambda col: col.map(lambda geneid: geneid.split('.')[1]))
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

pd.DataFrame([node for cond, node in zip(select, nodes) if cond]).to_csv('ensp.txt', index=False, header=None)
sp.save_npz('adj.npz', adj.tocsr()[select].tocsc()[:, select].tocoo())
```

## PPI figure

In addition, we draw the protein-protein interaction network of *E. coli* to intuitively show the
relation between the curvature and the network topology.
To reproduce the figure, run the following Python script:

```python
from collections import defaultdict

import numpy as np
import scipy.sparse as sp
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from pwn.utils import compute_curvature


edges = pd.read_csv(
    'https://stringdb-static.org/download/protein.physical.links.v11.0/511145.protein.physical.links.v11.0.txt.gz', sep=' '
).query('combined_score >= 700')
nodes = sorted(set(edges['protein1'].to_list()))
node2idx = {n: i for i, n in enumerate(nodes)}

adj = sp.dok_matrix((len(nodes), len(nodes)), bool)
for _, src, dst, _ in edges.itertuples():
    adj[node2idx[src], node2idx[dst]] = adj[node2idx[dst], node2idx[src]] = True

conn = sp.csgraph.connected_components(adj)[1]
n_conn = defaultdict(int)
for c in conn:
    n_conn[c] += 1

select = conn == max(n_conn, key=lambda x: n_conn[x])
nodes = ['.'.join(node.split('.')[1:]) for sel, node in zip(select, nodes) if sel]
adj = adj.tocsr()[select].tocsc()[:, select].tocoo()

G = nx.from_scipy_sparse_matrix(compute_curvature(adj), edge_attribute='curvature')
curv_mean = np.array([np.mean([e['curvature'] for e in G[i].values()]) for i in G])
curv = np.array([e[2]['curvature'] for e in G.edges.data()])

mm = 0.1 / 2.54
plt.figure(figsize=(0.5 * 170 * mm, 0.5 * 113 * mm))
plt.tight_layout(pad=0)
nx.draw(
    G,
    pos=nx.spring_layout(G, seed=31415),
    node_size=1.0,
    width=0.2,
    node_color=list(map(plt.cm.coolwarm, ((adj.sum(0).A1 >= 10) * np.arctan(0.5 * curv_mean) + 0.5 * np.pi) / np.pi)),
    edge_color=list(map(plt.cm.coolwarm, (np.arctan(0.5 * curv) + 0.5 * np.pi) / np.pi)),
)
plt.savefig('ppi.pdf', dpi=300)
```

## License

The data is originally provided by the STRING Consortium and redistributed under a
[Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).
