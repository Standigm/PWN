# Copyright 2022 Standigm Inc.. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

"""Official implementation of PWN and it's variants, with re-implemented other well-known methods."""

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from .utils import build_laplacian, compute_curvature, normalize, sigmoid, solve_linsys


def gdc(
    adj: sp.spmatrix, shift: float = 0.0, normalization: str = 'rw', restartprob: float = 0.12, topk: int = 128
) -> sp.spmatrix:
    r"""GDC: graph diffusion convolution (NeurIPS 2019).

    Note:
        Unlike other methods that calibrates gene scores, GDC generates new weighted adjacency matrix.

    Args:
        adj: Adjacency matrix :math:`\{0,1\}^{n\times n}`.
        shift: Weight of additional self-loop.
        normalization: Types of normalization; random walk normalize (``rw``) and symmetric normalize (``sym``).
        restartprob: Restart probability.
        topk: Number of edges to keep in sparsification.

    Returns:
        Weighted adjacency matrix :math:`\mathbb{R}_+^{n\times n}`.
    """
    wadj = normalize(adj.astype(float) + shift * sp.diags(np.ones(adj.shape[0])), normalization)
    S = restartprob * (sp.diags(np.ones(wadj.shape[0])) + (restartprob - 1) * wadj)
    S = np.linalg.inv(S.todense())
    index = np.argpartition(-S, topk, axis=0)[:topk, :]
    index = [(i, j) for j in range(S.shape[1]) for i in index[:, j].A1]
    return sp.coo_matrix(([S[i, j] for i, j in index], zip(*index)), shape=S.shape)


def rwr(adj: sp.spmatrix, x: npt.NDArray[np.float_], restartprob: float = 0.3) -> npt.NDArray[np.float_]:
    r"""Random walk with restart.

    Args:
        adj: Adjacency matrix :math:`\mathbb{R}_+^{n\times n}`; possibly weighted.
        x: Original gene score :math:`\mathbb{R}_+^n`.
        restartprob: Restart probability.

    Returns:
        Final gene score :math:`\mathbb{R}^n`.
    """
    return restartprob * solve_linsys(sp.eye(adj.shape[1]) - (1 - restartprob) * normalize(adj, 'rw'), x)


def mND(
    adj: sp.spmatrix, x: npt.NDArray[np.float_], restartprob: float = 0.3, numnbd: int = 3, eps: float = 1e-15
) -> npt.NDArray[np.float_]:
    r"""mND: (Bioinformatics 2020).

    Args:
        adj: Adjacency matrix :math:`\{0,1\}^{n\times n}`.
        x: Original gene scores :math:`\mathbb{R}^{n\times k}`
        restartprob: Restart probability used in score diffusion.
        numnbd: Number of neighbors used in score integration.
        eps: Small value to avoid numerical errors.

    Returns:
        Final gene score :math:`\mathbb{R}^n`.
    """
    transition = normalize(adj, 'sym')
    xnew = []
    for i in range(x.shape[1]):
        tmp = solve_linsys(sp.eye(transition.shape[0]) - (1 - restartprob) * transition, x[:, i])
        xnew.append(tmp / np.abs(tmp).max())
    xnew = np.stack(xnew, axis=-1)

    transition = transition.tocsr()
    T = np.zeros_like(xnew)
    for i in range(T.shape[0]):
        tmp = xnew[transition[i].indices]
        if transition[i].nnz > numnbd:
            tmp = -np.partition(-tmp, numnbd, axis=0)[:numnbd]
        T[i] = tmp.sum(axis=0)
    return T.sum(axis=-1) * xnew.sum(axis=-1) / (eps + np.minimum(numnbd, transition.astype(bool).sum(axis=-1).A1))


def uKIN(
    adj: sp.spmatrix,
    x: npt.NDArray[np.float_],
    prior: npt.NDArray[np.float_],
    shift: float = 1.0,
    restartprob: float = 0.3,
    eps: float = 1e-15,
) -> npt.NDArray[np.float_]:
    r"""uKIN: using knowledge in networks (RECOMB 2020).

    Args:
        adj: Adjacency matrix :math:`\{0,1\}^{n\times n}`.
        x: Original gene score :math:`\mathbb{R}_+^n`.
        prior: Indicator vector `\{0,1\}^n` that represents prior knowledge set of genes.
        shift: Amount of shift used in building the shifted Laplacian.
        restartprob: Restart probability.
        eps: Small value to avoid numerical errors.

    Returns:
        Final gene score :math:`\mathbb{R}^n`.
    """
    wadj = adj.astype(float, copy=True).tocsr(copy=False)
    prior = solve_linsys(build_laplacian(adj, shift), prior)
    for i in range(wadj.shape[0]):
        wadj.data[wadj.indptr[i] : wadj.indptr[i + 1]] *= prior[i] + eps
    return rwr(wadj, x, restartprob)


def pwn(
    adj: sp.spmatrix,
    x: npt.NDArray[np.float_],
    prior: npt.NDArray[np.float_],
    beta: float = 0.5,
    gamma: float = 0.5,
    restartprob: float = 0.3,
    precomputed: bool = True,
    eps: float = 1e-15,
) -> npt.NDArray[np.float_]:
    r"""PWN: prioritization with a warped network.

    Args:
        adj: Adjacency matrix :math:`\{0,1\}^{n\times n}`.
        x: Original gene score :math:`\mathbb{R}^n`
        prior: Indicator vector `\{0,1\}^n` that represents prior knowledge set of genes.
        beta: Amount of curvature shifting.
        gamma: Restart probability used in prior smoothing.
        restartprob: Restart probability used in gene score propagation.
        precomputed: Indicates whether curvatures are already computed or not.
        eps: Small value to avoid numerical errors.

    Returns:
        Final gene score :math:`\mathbb{R}^n`.
    """
    wadj = adj.astype(float, copy=True) if precomputed else compute_curvature(adj, inplace=False)
    wadj.data = sigmoid(beta * (wadj.data - wadj.data.mean()) / wadj.data.std())

    prior = rwr(wadj, prior / prior.sum(), restartprob=gamma)
    wadj = wadj.tocsr(copy=False)
    for i in range(wadj.shape[0]):
        wadj.data[wadj.indptr[i] : wadj.indptr[i + 1]] *= prior[i] + eps

    return rwr(wadj, x, restartprob=restartprob)
