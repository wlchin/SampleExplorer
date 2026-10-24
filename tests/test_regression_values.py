"""
Regression values for the numerical core: ORA (sample_explorer.ora), transcriptome
and semantic similarity (scikit-learn, sentence-transformers) and the search ranking.

The ORA values are decoupler 1.6's, which the package reproduces exactly. They were
captured on 2026-10-24; the other values with the uv.lock environment (pandas 3.1.0,
numpy 2.5.4, anndata 0.13.4, sentence-transformers 6.2.0, torch 2.14.1,
scikit-learn 1.9.1). A failure after a dependency upgrade means search results
changed; decide whether the change is acceptable before updating the expected values.

decoupler 2.3 was briefly used and differed from 1.6: its test is two-sided, and it
breaks ties between genes at the top-5% cutoff differently, which changed the p-value
of 13 of the first 500 fixture samples (15 of all 2,000). The semantic and
transcriptome similarity values did not change.
"""
import numpy as np
import pandas as pd
import pytest

from tests.common import TEST_GENESET, TEXT_QUERY


@pytest.fixture(scope="module")
def ora_first_500(query_db):
    te = query_db.transcriptome_enrichment
    return te.run_decouplr_on_memmaped_adata_with_samples(TEST_GENESET, te.memmap_adata.obs.index[:500])


@pytest.fixture(scope="module")
def ora_all(query_db):
    # same batching as Query_DB.transcriptome_search
    te = query_db.transcriptome_enrichment
    samples = te.memmap_adata.obs.index
    return pd.concat([te.run_decouplr_on_memmaped_adata_with_samples(TEST_GENESET, samples[i:i + 500])
                      for i in range(0, len(samples), 500)])


def test_ora_pvalue_distribution(ora_first_500):
    values, counts = np.unique(ora_first_500["pvals"], return_counts=True)
    np.testing.assert_allclose(values, [3.7538531891624514e-05, 0.002258355904199917, 0.07270375322600171, 1.0],
                               rtol=1e-6)
    assert counts.tolist() == [3, 28, 156, 313]


def test_ora_values_for_selected_samples(ora_first_500):
    # pvals, estimate (-log10 pvals), p_corrected (BH within this 500-sample batch)
    expected = pd.DataFrame(
        [[3.7538531891624514e-05, 4.425522716366618, 0.006256421981937418],
         [0.002258355904199917, 2.6462076146475404, 0.03642509522903092],
         [0.07270375322600171, 1.1384431687383787, 0.19439506210160887],
         [1.0, 0.0, 1.0]],
        index=["GSM1093230", "GSM1020215", "GSM1000981", "GSM1000984"],
        columns=["pvals", "estimate", "p_corrected"],
    )
    actual = ora_first_500.loc[expected.index, expected.columns]
    np.testing.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=1e-6, atol=1e-12)


def test_transcriptome_search_ranking(query_db, ora_all):
    res_df, series = query_db.transcriptome_search(TEST_GENESET, nsamples=10)
    assert len(res_df) == 10
    # Three samples score clearly highest. The other seven are drawn from a large
    # group tied at the next score, so which of those are returned depends on sort
    # order; only their score is pinned here.
    assert {"GSM1093230", "GSM1098362", "GSM1098366"} <= set(res_df.index)
    scores = sorted(ora_all.loc[res_df.index, "estimate"], reverse=True)
    np.testing.assert_allclose(scores, [4.425522716366618] * 3 + [2.6462076146475404] * 7, rtol=1e-6)
    assert series.tolist() == res_df["gse_id"].tolist()


def test_semantic_search_scores(query_db):
    df, _ = query_db.semantic_search(TEXT_QUERY, k=5)
    assert df["gse_id"].tolist() == ["GSE45157", "GSE46908", "GSE44615,GSE44616", "GSE44568,GSE44616", "GSE42951"]
    # query_rag reports cosine distance (lower is closer) in "similarity_score"
    np.testing.assert_allclose(df["similarity_score"],
                               [0.2946821451187134, 0.3558880090713501, 0.3953213095664978,
                                0.4235578775405884, 0.5480091571807861], rtol=1e-4)


def test_semantic_expansion_scores(query_db):
    df = query_db.rag_embedding.get_closest_semantic_studies("GSE45157", k=5)
    assert df["gse_id"].tolist() == ["GSE45157", "GSE45159", "GSE46908", "GSE44568,GSE44616", "GSE44615,GSE44616"]
    np.testing.assert_allclose(df["similarity_score"],
                               [1.000000238418579, 0.6481589078903198, 0.6433640718460083,
                                0.6253336668014526, 0.6253006458282471], rtol=1e-5)


@pytest.mark.parametrize("series_id, expected_ids, expected_scores", [
    ("GSE43801",
     ["GSE43801", "GSE43770", "GSE45326", "GSE40819,GSE40820", "GSE40819,GSE40820"],
     [0.9999998807907104, 0.6293402314186096, 0.5972867012023926, 0.5962581634521484, 0.5959545373916626]),
    ("GSE46323",
     ["GSE46323", "GSE46323", "GSE43143", "GSE47774,GSE47792", "GSE47774,GSE47792"],
     [0.9720855951309204, 0.9720855951309204, 0.896048903465271, 0.8659668564796448, 0.8656827211380005]),
])
def test_transcriptome_expansion_scores(query_db, series_id, expected_ids, expected_scores):
    df = query_db.transcriptome_embedding.get_closest_transcriptional_studies(series_id, k=5)
    assert df["gse_id"].tolist() == expected_ids
    np.testing.assert_allclose(df["similarity_score"], expected_scores, rtol=1e-5)


def test_sentence_embedding_similarity(query_db):
    sim = query_db.rag_embedding.get_averages_between_queries(
        "lincRNA regulation", ["long non-coding RNA", "cardiac myocytes"])
    np.testing.assert_allclose(sim, 0.1331131, rtol=1e-4)
