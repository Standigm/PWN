# Copyright 2022 Standigm Inc.. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

"""Collection of utility functions."""

import random

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from scipy.sparse.linalg import gcrotmk


def set_seed(seed: int):
    """Set seed."""
    random.seed(seed)
    np.random.seed(seed)


def sigmoid(x: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
    """Numerically stable sigmoid function."""
    return np.exp(x / 2) / (np.exp(x / 2) + np.exp(-x / 2))


def solve_linsys(A: sp.spmatrix, b: npt.NDArray[np.float_], tol: float = 1e-05) -> npt.NDArray[np.float_]:
    r"""Find the solution of :math:`Ax = b`.

    Args:
        A: :math:`\mathbb{R}^{n\times n}`.
        b: :math:`\mathbb{R}^n`.
        tol: Tolerances for convergence.

    Returns:
        :math:`\mathbb{R}^n`.
    """
    x, status_code = gcrotmk(A, b, tol=tol, atol=tol)
    assert status_code == 0, 'Not converged.'
    return x


def normalize(adj: sp.spmatrix, method: str = 'rw', eps: float = 1e-15, inplace: bool = False) -> sp.spmatrix:
    r"""Normalize the given adjacency matrix.

    Currently, two methods are available; the random walk normalization

    .. math::
       \bar{A} = AD^{-1},

    and the symmetric normalization

    .. math::
       \bar{A} = D^{-1/2}AD^{-1/2},

    where :math:`A` is a adjacency matrix and :math:`D` is the diagonal node degree matrix.

    Args:
        adj: Adjacency matrix :math:`\mathbb{R}_+^{n\times n}`; possibly weighted.
        method: Types of normalization; random walk normalize (``rw``) and symmetric normalize (``sym``).
        eps: Small value to avoid division by zero.
        inplace: Change original adjacency matrix or not.

    Returns:
        :math:`\mathbb{R}^{n\times n}`.
    """
    assert np.all(adj.data >= 0), 'Invalid adjacency matrix'
    D = adj.sum(axis=0, dtype=float).A1
    if method == 'rw':
        transition = adj.tocsc(copy=not inplace).astype(float)
        for i in range(transition.shape[1]):
            transition.data[transition.indptr[i] : transition.indptr[i + 1]] /= D[i] + eps
    elif method == 'sym':
        transition = adj.tocsr(copy=not inplace).astype(float)
        for i in range(transition.shape[0]):
            transition.data[transition.indptr[i] : transition.indptr[i + 1]] /= np.sqrt(D[i]) + eps
        transition = transition.tocsc(copy=False)
        for i in range(transition.shape[1]):
            transition.data[transition.indptr[i] : transition.indptr[i + 1]] /= np.sqrt(D[i]) + eps
    else:
        raise NotImplementedError('Possible types of normalization: rw (random walk), sym (symmetric).')
    return transition


def build_laplacian(adj: sp.spmatrix, delta: float = 0.0) -> sp.spmatrix:
    r"""Build the shifted Laplacian of given adjacency matrix.

    The graph Laplacian :math:`L` shifted by :math:`\delta` is

    .. math::
       L = D + \delta I - A,

    where :math:`A` is a adjacency matrix and :math:`D` is the diagonal node degree matrix.

    Args:
        adj: Adjacency matrix :math:`\mathbb{R}_+^{n\times n}`; possibly weighted.
        delta: Amount of shift.

    Returns:
        :math:`\mathbb{R}^{n\times n}`.
    """
    assert np.all(adj.data >= 0), 'Invalid adjacency matrix'
    assert delta >= 0.0, 'The amount of shift should be non-negative.'
    return sp.diags(adj.sum(axis=0, dtype=float).A1 + delta) - adj


def compute_curvature(adj: sp.spmatrix, inplace: bool = False) -> sp.spmatrix:
    r"""Compute the Forman-Ricci curvatures of every edges in network.

    Args:
        adj: Adjacency matrix :math:`\{0,1\}^{n\times n}`.
        inplace: Change original adjacency matrix or not.

    Returns:
        :math:`\mathbb{R}^{n\times n}`.
    """
    result = adj.tocoo(copy=not inplace).astype(int)
    edges = adj.tocsc(copy=True)
    for i, (v, w) in enumerate(zip(result.row, result.col)):
        v = edges.indices[edges.indptr[v] : edges.indptr[v + 1]]
        w = edges.indices[edges.indptr[w] : edges.indptr[w + 1]]
        result.data[i] = 4 - v.size - w.size + 3 * np.intersect1d(v, w, True).size
    return result.astype(float)
