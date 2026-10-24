"""
RNASeqAnalysis against a miniature ARCHS4 file (see tests/common.py), so the
ARCHS4 reader and the ORA run for real.

Expected values were captured on 2026-10-24 with archs4py 0.2.12 / decoupler 1.6.0,
which the package reproduces (sample_explorer.ora). decoupler 2.3 was briefly used
and differed (two-sided test, tie rule).
"""
import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control

from sample_explorer.rnaseq_analysis import RNASeqAnalysis
from tests.common import TEST_GENESET, read_archs4_counts

SAMPLES = ["GSM1000981", "GSM1000984", "GSM1002540", "GSM1009637", "GSM1020213", "GSM1093230", "GSM1158460"]


@pytest.fixture(scope="module")
def rnaseq(archs4_h5):
    return RNASeqAnalysis(archs4_h5)


@pytest.fixture(scope="module")
def enrichment(rnaseq):
    return rnaseq.perform_enrichment_on_samples(SAMPLES, TEST_GENESET)


@pytest.fixture(scope="module")
def enrichment_batched(rnaseq):
    return rnaseq.perform_enrichment_on_samples_batched(SAMPLES, TEST_GENESET, batch_size=4)


def test_create_anndata_from_samples(rnaseq, archs4_h5):
    adata = rnaseq.create_anndata_from_samples(SAMPLES)
    assert sorted(adata.obs_names) == sorted(SAMPLES)
    assert adata.n_vars == 122
    # metadata rows must line up with the expression rows they describe
    for sample_id in adata.obs_names:
        np.testing.assert_array_equal(adata[sample_id].X.ravel(), read_archs4_counts(archs4_h5, sample_id))


def test_create_anndata_from_series(rnaseq):
    adata = rnaseq.create_anndata_from_series("GSE29282")
    assert adata.obs_names.tolist() == ["GSM1000981", "GSM1000984"]
    assert (adata.obs["series_id"] == "GSE29282").all()


def test_list_to_dc_geneset(rnaseq):
    df = rnaseq.list_to_dc_geneset(TEST_GENESET)
    assert df.shape == (6, 2)
    assert df["genesymbol"].tolist() == TEST_GENESET
    assert set(df["geneset"]) == {"enrichment_score"}


def test_enrichment_columns(enrichment):
    assert enrichment.columns[:3].tolist() == ["enrichment_score", "pvals", "p_corrected"]
    assert {"geo_accession", "series_id", "extract_protocol_ch1"} <= set(enrichment.columns)
    assert sorted(enrichment.index) == sorted(SAMPLES)


def test_enrichment_values(enrichment):
    # enrichment_score is -log10(pvals); p_corrected is BH over these 7 samples
    expected = pd.DataFrame(
        [[2.6781064559859353, 0.002098425446207324, 0.002937795624690253],
         [2.6781064559859353, 0.002098425446207324, 0.002937795624690253],
         [0.0, 1.0, 1.0],
         [2.6781064559859353, 0.002098425446207324, 0.002937795624690253],
         [5.802987316603255, 1.5740288326778528e-06, 5.509100914372485e-06],
         [9.279971032072446, 5.248424666616609e-10, 3.6738972666316266e-09],
         [0.0, 1.0, 1.0]],
        index=SAMPLES,
        columns=["enrichment_score", "pvals", "p_corrected"],
    )
    actual = enrichment.loc[SAMPLES, expected.columns]
    np.testing.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=1e-6, atol=1e-12)


def test_single_batch_matches_unbatched(rnaseq, enrichment):
    batched = rnaseq.perform_enrichment_on_samples_batched(SAMPLES, TEST_GENESET)
    pd.testing.assert_frame_equal(batched.loc[SAMPLES], enrichment.loc[SAMPLES])


def test_batched_covers_all_samples(enrichment, enrichment_batched):
    assert sorted(enrichment_batched.index) == sorted(SAMPLES)
    assert enrichment_batched.columns.tolist() == enrichment.columns.tolist()


@pytest.mark.xfail(raises=AssertionError, reason="decoupler drops genes that are zero across the whole batch, so n_up (5% of genes) "
                          "and tie-breaking depend on which samples share a batch")
def test_scores_do_not_depend_on_batch_composition(enrichment, enrichment_batched):
    cols = ["enrichment_score", "pvals"]
    np.testing.assert_allclose(enrichment_batched.loc[SAMPLES, cols].to_numpy(), enrichment.loc[SAMPLES, cols].to_numpy())


@pytest.mark.xfail(raises=AssertionError, reason="BH correction is applied within each batch instead of across all samples")
def test_batched_fdr_is_corrected_across_all_samples(enrichment_batched):
    expected = false_discovery_control(enrichment_batched["pvals"], method="bh")
    np.testing.assert_allclose(enrichment_batched["p_corrected"], expected)


def test_perform_enrichment_on_series(rnaseq):
    res = rnaseq.perform_enrichment_on_series(["GSE29282", "GSE43801"], TEST_GENESET)
    assert set(res["series_id"]) == {"GSE29282", "GSE43801"}
    assert res.index.is_unique
    assert len(res) == 3
