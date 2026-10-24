"""
Checks the over-representation analysis (ORA) against an independent computation
with scipy, so that a result which is stable but wrong is caught (the regression
values only check that results do not change).

For each sample, the top 5% of expressed genes is tested for overlap with the gene
set using the hypergeometric upper tail (one-sided Fisher's exact test) against a
background of 20,000 genes, after excluding genes that are zero in every sample of
the batch. Which of several tied genes lands inside the top 5% is decided by
decoupler's internal shuffle, so only samples whose top 5% is unambiguous are
compared.
"""
import numpy as np
import pytest
from scipy.stats import false_discovery_control, hypergeom

from sample_explorer.rnaseq_analysis import RNASeqAnalysis
from tests.common import TEST_GENESET, read_archs4_matrix


def reference_ora(matrix, genes, geneset, n_background=20000):
    """Reference p-values, and a mask of samples with no tie at the top-5% cutoff."""
    matrix = np.asarray(matrix, dtype=np.float64)
    expressed = (matrix != 0).any(axis=0)
    matrix, genes = matrix[:, expressed], np.asarray(genes)[expressed]
    n_up = max(int(np.ceil(0.05 * len(genes))), 2)
    in_set = np.isin(genes, geneset)
    order = np.argsort(-matrix, axis=1, kind="stable")
    overlap = in_set[order[:, :n_up]].sum(axis=1)
    pvals = hypergeom.sf(overlap - 1, n_background, in_set.sum(), n_up)
    ranked = np.take_along_axis(matrix, order, axis=1)
    unambiguous = ranked[:, n_up - 1] > ranked[:, n_up]
    return pvals, unambiguous


@pytest.fixture(scope="module")
def transcriptome_batch(query_db):
    te = query_db.transcriptome_enrichment
    samples = te.memmap_adata.obs.index[:500]
    result = te.run_decouplr_on_memmaped_adata_with_samples(TEST_GENESET, samples)
    # same preprocessing as run_decouplr_on_memmaped_adata_with_samples (float16 overflow becomes 0)
    matrix = np.nan_to_num(te.memmap_adata.X[:500], nan=0.0, posinf=0.0, neginf=0.0)
    return result.loc[samples], matrix, te.memmap_adata.var_names


@pytest.fixture(scope="module")
def rnaseq_batch(archs4_h5):
    genes, samples, matrix = read_archs4_matrix(archs4_h5, n_samples=300)
    result = RNASeqAnalysis(archs4_h5).perform_enrichment_on_samples(samples, TEST_GENESET)
    return result.loc[samples], matrix, genes


def test_transcriptome_ora_pvalues_match_reference(transcriptome_batch):
    result, matrix, genes = transcriptome_batch
    expected, unambiguous = reference_ora(matrix, genes, TEST_GENESET)
    # enough comparable samples, including enriched ones, for the check to mean something
    assert unambiguous.sum() >= 200
    assert (unambiguous & (expected < 1)).sum() >= 50
    np.testing.assert_allclose(result["pvals"].to_numpy()[unambiguous], expected[unambiguous], rtol=1e-6)


def test_transcriptome_ora_estimate_and_fdr(transcriptome_batch):
    result, _, _ = transcriptome_batch
    np.testing.assert_allclose(result["estimate"], -np.log10(result["pvals"]))
    np.testing.assert_allclose(result["p_corrected"], false_discovery_control(result["pvals"], method="bh"))


def test_rnaseq_ora_pvalues_match_reference(rnaseq_batch):
    result, matrix, genes = rnaseq_batch
    expected, unambiguous = reference_ora(matrix, genes, TEST_GENESET)
    assert unambiguous.sum() >= 150
    assert (unambiguous & (expected < 1)).sum() >= 50
    np.testing.assert_allclose(result["pvals"].to_numpy()[unambiguous], expected[unambiguous], rtol=1e-6)


def test_rnaseq_ora_estimate_and_fdr(rnaseq_batch):
    result, _, _ = rnaseq_batch
    np.testing.assert_allclose(result["enrichment_score"], -np.log10(result["pvals"]))
    np.testing.assert_allclose(result["p_corrected"], false_discovery_control(result["pvals"], method="bh"))
