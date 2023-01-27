# Copyright 2022 Standigm Inc.. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

"""Reproducible experiments for the paper."""

import os
import re
from argparse import ArgumentParser
from enum import Enum
from typing import Any, Optional, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse as sp
import scipy.stats as ss
from sklearn.metrics import average_precision_score
from tqdm import tqdm

import pwn.method as methods
from pwn.utils import build_laplacian, compute_curvature, set_seed, sigmoid, solve_linsys

BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
METHOD_REGEX = re.compile(r'^RWR \+ Prior\((uKIN|RWR)(, Curv)?\)( \+ Curv)?$')


class GeneID(Enum):
    """Available ID system for genes."""

    ENSG = 'ensg'
    ENSP = 'ensp'
    UniProt = 'uniprot'


class PPI(Enum):
    """Available sources of protein-protein interaction networks."""

    STRING = 'string'
    BioGRID = 'biogrid'
    IID = 'iid'

    def __str__(self):
        return self.value


class Result(Enum):
    """Result files."""

    TOPOLOGY = 'topo'
    SCORE = 'score'
    PERFORMANCE = 'perf'
    VARIANCE = 'var'

    def save(self, obj: pd.DataFrame, source: PPI):
        """Save to tsv."""
        obj.to_csv(f'{self.value}.{source}.tsv', index=False, sep='\t')


def compute_zscore(pvalue: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
    r"""Compute :math:`Z`-score from :math:`p`-value."""
    with np.errstate(invalid='ignore'):
        return np.nan_to_num(ss.norm.isf(pvalue), nan=0.0, posinf=ss.norm.isf(np.finfo(float).tiny), neginf=0.0)


def compute_nlogp(pvalue: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
    r"""Compute :math:`-\log(p)`."""
    with np.errstate(divide='ignore'):
        return np.nan_to_num(-np.log(pvalue), nan=0.0, posinf=-np.log(np.finfo(float).tiny))


def combine_pvalues(pvalues: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
    r"""Combine :math:`p`-values via Fisher method."""
    num_valids = np.isfinite(pvalues).sum(axis=-1)
    nlogp = compute_nlogp(pvalues).sum(axis=-1)
    return np.nan_to_num(ss.chi2.sf(2 * nlogp, 2 * num_valids), nan=0.0)


def load_ppi(source: PPI = PPI.STRING) -> tuple[GeneID, list[str], sp.spmatrix]:
    r"""Load protein-protein interaction network, with corresponding ENSP IDs.

    Args:
        source: Source of the network.

    Returns:
        A tuple containing:
        - ID system; ENSG or ENSP.
        - List of IDs
        - Protein-protein interaction network :math:`\{0,1\}^{m\times m}`.
    """
    if source == PPI.STRING:
        idtype = GeneID.ENSP
    elif source == PPI.BioGRID:
        idtype = GeneID.ENSG
    elif source == PPI.IID:
        idtype = GeneID.UniProt
    else:
        raise NotImplementedError
    nodeid = pd.read_csv(os.path.join(BASE_DIR, 'ppi', source.value, f'{idtype.value}.txt'), header=None)[0].tolist()
    adj = sp.load_npz(os.path.join(BASE_DIR, 'ppi', source.value, 'adj.npz'))
    return idtype, nodeid, adj


def load_converters(source: PPI, idtype: GeneID) -> Optional[dict[str, str]]:
    r"""Load gene ID converter.

    Args:
        source: Source of the network.
        idtype: Target gene ID.

    Returns:
        Gene ID converter.
    """
    if idtype == GeneID.ENSG:
        return None
    df = pd.read_csv(os.path.join(BASE_DIR, 'idmap', f'{source}.tsv'), sep='\t', header=None)
    id_converter = {}
    for _, ensg, identifier in df.dropna().drop_duplicates().itertuples():
        id_converter[ensg] = identifier
    return id_converter


def load_prior(nodeid: list[str], converter: Optional[dict[str, str]]) -> sp.spmatrix:
    r"""Load prior genes, alinged by given list of IDs.

    Args:
        nodeid: List of IDs.
        converter: ID converter.

    Returns:
        Prior genes indicator :math:`\{0,1\}^{k\times m}`.
    """
    id2idx = {id: idx for idx, id in enumerate(nodeid)}
    priors = sp.dok_matrix((1, len(nodeid)), dtype=bool)

    for _, prior in pd.read_csv(os.path.join(BASE_DIR, 'prior', 'ensg.txt'), header=None).itertuples():
        if converter is not None:
            if prior in converter and converter[prior] in id2idx:
                priors[0, id2idx[converter[prior]]] = True
        else:
            if prior in id2idx:
                priors[0, id2idx[prior]] = True
    return priors.tocsr()


def load_omics(nodeid: list[str], converter: Optional[dict[str, str]]) -> npt.NDArray[np.float_]:
    r"""Load :math:`p`-value from omics data, alinged by given list of ENSP IDs.

    Args:
        nodeid: List of IDs.
        converter: ID converter.

    Returns:
        :math:`p`-values :math:`\mathbb{R}^{k\times m\times d}`.
    """
    omics = pd.read_csv(os.path.join(BASE_DIR, 'omics', 'pvalue.csv'), index_col='ensg')
    if converter is not None:
        converter_df = pd.DataFrame(converter.items(), columns=['ensg', '']).set_index('ensg')
        omics = converter_df.join(omics, how='left').reset_index().set_index('').drop('ensg', axis=1)
    df = pd.DataFrame(index=nodeid).join(omics, how='left')
    return df[~df.index.duplicated()].loc[nodeid].to_numpy()[np.newaxis, ...]


def load_and_preproc(
    ppi_source: PPI, num_replicates: int = 1, train_ratios: list[float] = [0.2]
) -> tuple[
    list[str],
    dict[str, sp.spmatrix],
    dict[str, Union[dict[float, npt.NDArray[np.float_]], npt.NDArray[np.float_]]],
    dict[str, npt.NDArray[np.float_]],
]:
    r"""Load and preprocess raw data.

    Args:
        ppi_source: Source of PPI; ``string``, ``biogrid``, or ``iid``.
        num_replicates: Number of replicated experiments.
        train_ratios: List of ratio of train set.
    """
    idtype, nodeid, ppi = load_ppi(ppi_source)
    id_converter = load_converters(ppi_source, idtype)
    prior = load_prior(nodeid, id_converter)
    pvalue = load_omics(nodeid, id_converter)

    adj = {'unweighted': ppi.astype(float), 'curvature': compute_curvature(ppi)}

    test_ratio = 1.0 - max(train_ratios)
    priors = sp.vstack([prior] * num_replicates)
    test = np.zeros(priors.shape, dtype=bool)
    train = {rho: np.zeros_like(test) for rho in train_ratios}
    for i, prior in enumerate(priors.tocsr()):
        test_select = np.random.choice(prior.indices, size=np.floor(prior.nnz * test_ratio).astype(int), replace=False)
        test[i, test_select] = True
        for j, rho in enumerate(sorted(train_ratios, reverse=True)):
            if j == 0:
                train_select = np.setdiff1d(prior.indices, test_select)
            else:
                train_select = np.random.choice(
                    np.setdiff1d(prior.indices, test_select), size=np.floor(prior.nnz * rho).astype(int), replace=False
                )
            train[rho][i, train_select] = True
    y = {'all': priors.A, 'train': train, 'test': test}

    pvalue = np.concatenate([pvalue] * num_replicates, axis=0)
    pvalue_combined = combine_pvalues(pvalue)
    X = {
        'zscore': compute_zscore(pvalue),
        'zscore_combined': compute_zscore(pvalue_combined),
        'nlogp': compute_nlogp(pvalue),
        'nlogp_combined': compute_nlogp(pvalue_combined),
    }

    return nodeid, adj, y, X


def run_exp(
    adj: dict[Any, sp.spmatrix],
    X: dict[str, npt.NDArray[np.float_]],
    y_train: Optional[npt.NDArray[np.float_]],
    method: str,
    beta: float = 0.5,
    gamma: float = 0.5,
    shift: float = 1.0,
) -> npt.NDArray[np.float_]:
    """Get new gene scores."""
    match = METHOD_REGEX.match(method)
    if match is not None:
        assert y_train is not None
        adj_prior = adj[0.0] if match.group(2) is None else adj[beta]
        if match.group(1) == 'uKIN':
            prior = solve_linsys(build_laplacian(adj_prior, shift), y_train / y_train.sum())
        else:
            prior = methods.rwr(adj_prior, y_train / y_train.sum(), restartprob=gamma)
        adj_score = (adj[0.0] if match.group(3) is None else adj[beta]).tocsr(copy=True)
        for i in range(adj_score.shape[0]):
            adj_score.data[adj_score.indptr[i] : adj_score.indptr[i + 1]] *= prior[i] + 1e-15
        result = methods.rwr(adj_score, X['zscore_combined'])
    else:
        if method == 'Naive':
            result = X['zscore_combined']
        elif method == 'RWR':
            result = methods.rwr(adj['unweighted'], X['zscore_combined'])
        elif method == 'RWR w/ GDC':
            result = methods.rwr(adj['GDC'], X['zscore_combined'])
        elif method == 'RWR + Curv':
            result = methods.rwr(adj[beta], X['zscore_combined'])
        elif method == 'mND':
            result = methods.mND(adj['unweighted'], X['zscore'])
        else:
            raise NotImplementedError(method)
    return result


def measure_perf(
    y_train: npt.NDArray[np.float_],
    y_test: npt.NDArray[np.float_],
    y_pred: npt.NDArray[np.float_],
    ks: list[int] = [100, 200],
) -> dict[str, float]:
    """Measure several performance metrics."""
    filtered_y_test, filtered_y_pred = y_test[np.logical_not(y_train)], y_pred[np.logical_not(y_train)]
    hits_at_k = np.cumsum(filtered_y_test[np.argsort(-filtered_y_pred)][: max(ks)])
    perfs = {f'Prec@{k}': hits_at_k[k - 1] / k for k in ks}
    perfs['AveP'] = average_precision_score(filtered_y_test, filtered_y_pred)
    return perfs


def collect_experiment_result(args, adj, y, X, ks=[100, 200]):
    """Run a single experiment and return result tuples."""
    idx, fn, rho, beta, gamma = args
    kwarg = {}
    if np.isfinite(beta):
        kwarg['beta'] = beta
    if np.isfinite(gamma):
        kwarg['gamma'] = gamma
    if np.isfinite(rho):
        y_train, y_test = y['train'][rho][idx], y['test'][idx]
    else:
        y_train, y_test = np.zeros_like(y['all'][idx]), y['all'][idx]
    return [
        (idx, rho, beta, gamma, fn, k, v)
        for k, v in measure_perf(y_train, y_test, run_exp(adj, X, y_train, fn, **kwarg), ks=ks).items()
    ]


def prior_smoothing(
    adj: dict[str, sp.spmatrix],
    y_train: npt.NDArray[np.float_],
    method: str,
    beta: float = 0.5,
    gamma: float = 0.5,
    shift: float = 1.0,
) -> npt.NDArray[np.float_]:
    """Smoothed prior knowledge."""
    if method == 'uKIN':
        return solve_linsys(build_laplacian(adj['unweighted'], shift), y_train / y_train.sum())
    elif method == 'RWR':
        wadj = adj['curvature'].astype(float, copy=True)
        wadj.data = sigmoid(beta * (wadj.data - wadj.data.mean()) / wadj.data.std())
        return methods.rwr(wadj, y_train / y_train.sum(), restartprob=gamma)
    else:
        raise NotImplementedError


def main(seed: int, ppi: PPI):
    """Entrypoint."""
    num_replicates = 30
    rho0, rhos = 0.2, np.linspace(0.20, 0.05, 4).tolist()
    beta0, betas = 0.5, np.linspace(-1.0, 1.0, 5).tolist()
    gamma0, gammas = 0.5, np.linspace(0.1, 0.9, 5).tolist()

    set_seed(seed)
    nodeid, adj, y, X = load_and_preproc(ppi, num_replicates=num_replicates, train_ratios=rhos)
    adj['GDC'] = methods.gdc(adj['unweighted'])
    for beta in betas:
        if beta != 0.0:
            adj[beta] = adj['curvature'].copy()
            adj[beta].data = sigmoid(beta * (adj[beta].data - adj[beta].data.mean()) / adj[beta].data.std())
        else:
            adj[beta] = adj['unweighted'].astype(float, copy=True)
    X = {k: v[0] for k, v in X.items()}

    df = pd.DataFrame(
        {
            'id': nodeid,
            'prior': y['all'][0],
            'degree': adj['unweighted'].sum(axis=0).A1,
            'curv_mean': [r.data.mean() for r in adj['curvature'].tocsr()],
            'curv_std': [r.data.std() for r in adj['curvature'].tocsr()],
        }
    )
    Result.TOPOLOGY.save(df, ppi)

    exp_list = [
        # exp 1
        *[
            (i, fn, rho0, -np.inf, -np.inf)
            for i in range(num_replicates)
            for fn in ['Naive', 'RWR w/ GDC', 'mND', 'RWR']
        ],
        *[
            (i, fn, rho0, beta0, -np.inf)
            for i in range(num_replicates)
            for fn in [
                'RWR + Curv',
                'RWR + Prior(uKIN, Curv)',
                'RWR + Prior(uKIN) + Curv',
                'RWR + Prior(uKIN, Curv) + Curv',
            ]
        ],
        *[(i, fn, rho0, -np.inf, gamma0) for i in range(num_replicates) for fn in ['RWR + Prior(RWR)']],
        *[
            (i, fn, rho0, beta0, gamma0)
            for i in range(num_replicates)
            for fn in ['RWR + Prior(RWR, Curv)', 'RWR + Prior(RWR) + Curv']
        ],
        # exp 2
        *[(0, fn, -np.inf, -np.inf, -np.inf) for fn in ['Naive', 'RWR w/ GDC', 'mND']],
        *[(0, 'RWR + Curv', -np.inf, beta, -np.inf) for beta in betas],
        # exp 3
        *[
            (i, 'RWR + Prior(RWR, Curv) + Curv', rho0, beta, gamma)
            for i in range(num_replicates)
            for beta in betas
            for gamma in gammas
        ],
        # exp 4
        *[(i, 'RWR + Prior(uKIN)', rho, -np.inf, -np.inf) for i in range(num_replicates) for rho in rhos],
        *[
            (i, 'RWR + Prior(RWR, Curv) + Curv', rho, beta0, gamma0)
            for i in range(num_replicates)
            for rho in rhos
            if rho != rho0
        ],
    ]
    df = pd.DataFrame(
        sum(map(lambda args: collect_experiment_result(args, adj, y, X), tqdm(exp_list)), []),
        columns=['idx', 'rho', 'beta', 'gamma', 'method', 'metric', 'perf'],
    )
    Result.PERFORMANCE.save(df, ppi)

    scores = {
        'RWR': run_exp(adj, X, None, 'RWR'),
        'RWR + Curv': run_exp(adj, X, None, 'RWR + Curv', beta=beta0),
        'RWR + Prior(uKIN)': run_exp(adj, X, y['all'][0], 'RWR + Prior(uKIN)'),
        'RWR + Prior(uKIN, Curv)': run_exp(adj, X, y['all'][0], 'RWR + Prior(uKIN, Curv)', beta=beta0),
        'RWR + Prior(uKIN) + Curv': run_exp(adj, X, y['all'][0], 'RWR + Prior(uKIN) + Curv', beta=beta0),
        'RWR + Prior(uKIN, Curv) + Curv': run_exp(adj, X, y['all'][0], 'RWR + Prior(uKIN, Curv) + Curv', beta=beta0),
        'RWR + Prior(RWR)': run_exp(adj, X, y['all'][0], 'RWR + Prior(RWR)'),
        'RWR + Prior(RWR, Curv)': run_exp(adj, X, y['all'][0], 'RWR + Prior(RWR, Curv)', beta=beta0),
        'RWR + Prior(RWR) + Curv': run_exp(adj, X, y['all'][0], 'RWR + Prior(RWR) + Curv'),
        'RWR + Prior(RWR, Curv) + Curv': run_exp(adj, X, y['all'][0], 'RWR + Prior(RWR, Curv) + Curv', beta=beta0),
    }
    df = pd.DataFrame({'id': nodeid, **scores})
    Result.SCORE.save(df, ppi)

    variances = {
        method: np.stack(
            [
                prior_smoothing(adj, y['train'][rho0][idx], method, beta=beta0, gamma=gamma0)
                for idx in range(num_replicates)
            ]
        ).std(axis=0)
        for method in ['uKIN', 'RWR']
    }
    df = pd.DataFrame({'id': nodeid, **variances})
    Result.VARIANCE.save(df, ppi)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--seed', type=int, default=123)
    parser.add_argument('--ppi', type=PPI, choices=PPI, nargs='+', default=list(PPI))
    args = parser.parse_args()
    for ppi in args.ppi:
        main(args.seed, ppi)
