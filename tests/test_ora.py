import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom

from sample_explorer.ora import N_BACKGROUND, bh_adjust, run_ora
from tests.common import TEST_GENESET, TESTS_DIR, TRANSCRIPTOME_DB, read_archs4_matrix

DATA = TESTS_DIR / "data"


def _random_matrix(n_samples, n_genes, seed=0):
    """Positive continuous values: no gene is zero in every sample and there are no ties."""
    rng = np.random.default_rng(seed)
    matrix = rng.random((n_samples, n_genes)) + 0.01
    genes = [f"G{i}" for i in range(n_genes)]
    samples = [f"S{i}" for i in range(n_samples)]
    return matrix, samples, genes


def test_matches_decoupler16_on_transcriptome_fixture():
    expected = pd.read_csv(DATA / "decoupler16_transcriptome_pvals.csv", index_col="sample")["pvals"]
    adata = ad.read_h5ad(TRANSCRIPTOME_DB, backed="r")
    try:
        parts = []
        for i in range(0, adata.n_obs, 500):
            samples = adata.obs.index[i:i + 500]
            batch = adata[adata.obs.index.isin(samples)].to_memory()
            batch.var_names_make_unique()
            np.nan_to_num(batch.X, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
            pvals, _ = run_ora(batch.X, batch.obs.index, batch.var_names, {"query_gene_set": TEST_GENESET})
            parts.append(pvals["query_gene_set"])
    finally:
        adata.file.close()
    actual = pd.concat(parts)
    assert len(expected) == 2000
    assert list(actual.index) == list(expected.index)
    np.testing.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=1e-6)


def test_matches_decoupler16_on_mini_archs4(archs4_h5):
    expected = pd.read_csv(DATA / "decoupler16_rnaseq_pvals.csv", index_col="sample")["pvals"]
    genes, samples, matrix = read_archs4_matrix(archs4_h5, n_samples=300)
    pvals, _ = run_ora(matrix, samples, genes, {"enrichment_score": TEST_GENESET})
    assert len(expected) == 300
    assert list(pvals.index) == list(expected.index)
    np.testing.assert_allclose(pvals["enrichment_score"].to_numpy(), expected.to_numpy(), rtol=1e-6)


def test_independent_of_input_gene_order():
    matrix, samples, genes = _random_matrix(20, 300)
    gene_set = {"s": genes[:10] + genes[100:110]}
    p_original, _ = run_ora(matrix, samples, genes, gene_set)
    perm = np.random.default_rng(1).permutation(len(genes))
    p_shuffled, _ = run_ora(matrix[:, perm], samples, [genes[i] for i in perm], gene_set)
    pd.testing.assert_frame_equal(p_shuffled, p_original)


def test_all_zero_sample_has_p_one():
    matrix, samples, genes = _random_matrix(5, 100)
    matrix[2] = 0
    samples[2] = "zero_sample"
    pvals, _ = run_ora(matrix, samples, genes, {"s": genes[:10]})
    assert pvals.loc["zero_sample", "s"] == 1.0


def test_set_all_zero_in_batch_has_p_one():
    matrix, samples, genes = _random_matrix(5, 100)
    matrix[:, :10] = 0  # the set genes are in `genes` but zero in every sample
    pvals, _ = run_ora(matrix, samples, genes, {"s": genes[:10]})
    assert (pvals["s"] == 1.0).all()


def _run_single(gene_set):
    matrix, samples, genes = _random_matrix(8, 100, seed=3)
    pvals, scores = run_ora(matrix, samples, genes, {"s": gene_set})
    return pvals


def test_duplicate_set_genes_collapse():
    pd.testing.assert_frame_equal(_run_single(["G1", "G1", "G2"]), _run_single(["G1", "G2"]))


def test_absent_genes_ignored_and_fully_absent_set_raises():
    pd.testing.assert_frame_equal(
        _run_single(["G1", "G2", "NOT_A_GENE"]), _run_single(["G1", "G2"])
    )
    matrix, samples, genes = _random_matrix(3, 20)
    with pytest.raises(ValueError, match="empty_set"):
        run_ora(matrix, samples, genes, {"empty_set": ["NOT_A_GENE"]})


def test_duplicate_gene_names_rejected():
    matrix = np.ones((2, 3))
    with pytest.raises(ValueError, match="duplicate"):
        run_ora(matrix, ["x", "y"], ["A", "A", "B"], {"s": ["A"]})


def test_one_sided_for_depleted_samples():
    # 4000 genes, so n_up = 200; a set of 2000 genes has an expected overlap of 20.
    n_genes, n_set, n_up = 4000, 2000, 200
    rng = np.random.default_rng(5)
    matrix = rng.random((6, n_genes)) + 0.01  # strictly positive: no gene is dropped
    genes = [f"G{i}" for i in range(n_genes)]
    gene_set = genes[:n_set]
    # sample 0: top 200 contains no set gene; sample 1: top 200 is entirely set genes
    matrix[0, :n_set] = rng.random(n_set) * 0.5
    matrix[0, n_set:] = 1 + rng.random(n_genes - n_set)
    matrix[1, :n_set] = 1 + rng.random(n_set)
    matrix[1, n_set:] = rng.random(n_genes - n_set) * 0.5
    samples = ["depleted", "enriched", "a", "b", "c", "d"]
    pvals, _ = run_ora(matrix, samples, genes, {"s": gene_set})
    assert pvals.loc["depleted", "s"] == 1.0
    assert pvals.loc["enriched", "s"] < 1e-10
    in_set = np.isin(genes, gene_set)
    for i, sample in enumerate(samples):
        top = np.argsort(-matrix[i])[:n_up]
        a = in_set[top].sum()
        assert pvals.loc[sample, "s"] == pytest.approx(
            hypergeom.sf(a - 1, N_BACKGROUND, n_set, n_up), rel=1e-9
        )


def test_multiple_sets_are_raw_and_sorted():
    matrix, samples, genes = _random_matrix(10, 200, seed=7)
    set_a, set_b = genes[:8], genes[50:70]
    pvals, _ = run_ora(matrix, samples, genes, {"b": set_b, "a": set_a})
    assert list(pvals.columns) == ["a", "b"]
    assert pvals.dtypes.eq(np.float64).all()
    assert list(pvals.index) == samples
    only_a, _ = run_ora(matrix, samples, genes, {"a": set_a})
    only_b, _ = run_ora(matrix, samples, genes, {"b": set_b})
    pd.testing.assert_series_equal(pvals["a"], only_a["a"])
    pd.testing.assert_series_equal(pvals["b"], only_b["b"])


def test_scores_are_finite():
    # 40000 genes (n_up = 2000); the set is exactly the top 2000 genes, so p underflows to 0
    n_genes = 40000
    rng = np.random.default_rng(9)
    matrix = (rng.random((1, n_genes)) + 0.01)
    matrix[0, :2000] += 10
    genes = [f"G{i}" for i in range(n_genes)]
    pvals, scores = run_ora(matrix, ["s0"], genes, {"s": genes[:2000]})
    assert pvals.iloc[0, 0] == 0.0
    assert np.isfinite(scores.to_numpy()).all()
    assert scores.iloc[0, 0] == pytest.approx(-np.log10(np.finfo(float).tiny))


def test_bh_adjust():
    np.testing.assert_allclose(
        bh_adjust([0.01, 0.04, 0.03, 0.2]), [0.04, 0.16 / 3, 0.16 / 3, 0.2]
    )
