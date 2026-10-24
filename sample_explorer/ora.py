"""
Over-representation analysis (ORA) that reproduces decoupler 1.6 exactly.

For each sample the top 5% of genes by expression is tested for overlap with each
gene set using the upper tail of the hypergeometric distribution (a one-sided
Fisher's exact test) against a fixed background. Which of several tied genes falls
inside the top 5% follows decoupler 1.6's seeded permutation, so p-values match
its output, including on data with many ties (e.g. zero-inflated counts).
"""
from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, hypergeom

N_BACKGROUND = 20000
_CHUNK = 256
_TINY = np.finfo(np.float64).tiny


def _select_top(X, n_up):
    """Boolean mask of the n_up highest values per row; ties go to the later column."""
    t = np.partition(X, -n_up, axis=1)[:, [-n_up]]  # n_up-th largest per row
    above, tied = X > t, X == t
    slots = n_up - above.sum(axis=1, keepdims=True)
    from_right = np.cumsum(tied[:, ::-1], axis=1)[:, ::-1]
    return above | (tied & (from_right <= slots))


def run_ora(matrix, samples, genes, gene_sets: Mapping[str, Iterable[str]],
            n_background: int = N_BACKGROUND):
    """
    Over-representation analysis of each sample against each gene set.

    Args:
        matrix: samples x genes expression values (numpy, h5py or AnnData ``.X``).
        samples: sample identifiers, one per row of ``matrix``.
        genes: unique gene names, one per column of ``matrix``.
        gene_sets: mapping of gene-set name to its genes.
        n_background: number of genes in the background population.

    Returns:
        ``(pvals, scores)``: DataFrames indexed by ``samples`` with one float64
        column per gene set, sorted by name. ``pvals`` are raw (unadjusted) and
        ``scores`` are ``-log10(pvals)``, with p floored at the smallest normal float.
    """
    X = np.asarray(matrix, dtype=np.float32)  # decoupler 1.6 ranks float32 values
    genes = np.asarray(genes, dtype=str)
    if np.unique(genes).size != genes.size:
        raise ValueError("duplicate gene names in input")
    samples = list(samples)
    names = sorted(gene_sets)
    sets = {name: np.unique(np.asarray(list(gene_sets[name]), dtype=str)) for name in names}
    for name in names:
        if not np.isin(sets[name], genes).any():
            raise ValueError(f"gene set '{name}' has no genes in the input")

    keep = np.count_nonzero(X, axis=0) > 0  # drop genes that are zero in every sample
    X, genes = X[:, keep], genes[keep]
    order = np.argsort(genes)
    X, genes = X[:, order], genes[order]
    perm = np.arange(genes.size)
    np.random.default_rng(42).shuffle(perm)  # decoupler 1.6 tie-break permutation
    X, genes = X[:, perm], genes[perm]

    pvals = np.ones((len(samples), len(names)))
    if genes.size:
        n_up = int(np.ceil(0.05 * genes.size))
        in_sets = [np.isin(genes, sets[name]) for name in names]
        sizes = [int(m.sum()) for m in in_sets]
        for start in range(0, len(samples), _CHUNK):
            chunk = X[start:start + _CHUNK]
            selected = _select_top(chunk, n_up)
            nonzero = np.count_nonzero(chunk, axis=1) > 0
            for j, (in_set, k) in enumerate(zip(in_sets, sizes)):
                if k == 0:
                    continue
                a = (selected & in_set).sum(axis=1)
                p = hypergeom.sf(a - 1, n_background, k, n_up)
                pvals[start:start + _CHUNK, j] = np.where(nonzero, p, 1.0)

    scores = -np.log10(np.maximum(pvals, _TINY))
    index = pd.Index(samples)
    return (pd.DataFrame(pvals, index=index, columns=names),
            pd.DataFrame(scores, index=index, columns=names))


def bh_adjust(pvals) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    return false_discovery_control(pvals, method="bh")
