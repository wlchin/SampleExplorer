"""
End-to-end tests of Query_DB.search on the fixture stores, without and with an
ARCHS4 file (the miniature file built in tests/common.py).
"""
import logging

import pandas as pd
import pytest

from sample_explorer.archs4 import SAMPLE_FIELDS

text_query = "Studies that investigate severity of COVID-19 in adults"
# only IFNG, CCL2 and STAT1 are in the fixture stores
geneset_query = ["ACE2", "TMPRSS2", "IFITM3", "DPP4", "CTSL", "IL6", "TNF", "IFNG", "CCL2", "STAT1"]


@pytest.fixture(scope="module")
def geneset_search(query_db):
    return query_db.search(geneset=geneset_query, text_query=None, n_seed=20, n_expansion=1, perform_enrichment=False)


@pytest.fixture(scope="module")
def enriched_search(query_db_archs4):
    return query_db_archs4.search(geneset=geneset_query, text_query=None, n_seed=20, n_expansion=1,
                                  perform_enrichment=True)


@pytest.fixture(scope="module")
def metadata_search(query_db_archs4):
    return query_db_archs4.search(geneset=geneset_query, text_query=None, n_seed=20, n_expansion=1,
                                  perform_enrichment=False)


def test_search_working_seed_study_dimension(query_db):
    result = query_db.search(geneset=geneset_query, text_query=text_query, n_seed=20, n_expansion=0,
                             perform_enrichment=True)
    assert result.seed_studies.shape[0] == 20
    # without an ARCHS4 file enrichment is skipped
    assert result.samples is None


def test_search_working_should_be_none_if_n_expansion_none(query_db):
    result = query_db.search(geneset=geneset_query, text_query=text_query, n_seed=20, n_expansion=0,
                             perform_enrichment=False)
    assert result.expansion_studies is None


def test_text_only_search_is_semantic(query_db):
    result = query_db.search(geneset=None, text_query=text_query, n_seed=20, n_expansion=0, perform_enrichment=False)
    assert result.seed_studies.shape[0] == 20
    assert set(result.seed_studies["similarity_source"]) == {"semantic"}
    assert result.expansion_studies is None


def test_text_only_search_with_expansion(query_db):
    result = query_db.search(geneset=None, text_query=text_query, n_seed=20, n_expansion=1, perform_enrichment=False)
    assert result.expansion_studies.shape[0] == 20


def test_enrichment_without_geneset_raises(query_db):
    with pytest.raises(ValueError):
        query_db.search(geneset=None, text_query=text_query, n_seed=20, n_expansion=1, perform_enrichment=True)


def test_geneset_only_search_with_expansion(geneset_search):
    # transcriptome search returns 20 samples; each expands to k=1 study
    assert geneset_search.seed_studies.shape[0] == 20
    assert geneset_search.expansion_studies.shape[0] == 20
    assert geneset_search.samples is None


def test_seed_studies_are_part_of_expansion(geneset_search):
    seed_ids = geneset_search.seed_studies["gse_id"]
    assert seed_ids.isin(geneset_search.expansion_studies["gse_id"]).all()


def test_enrichment_adds_statistics(enriched_search):
    samples = enriched_search.samples
    assert samples.shape[0] > 0
    assert {"enrichment_score", "pvals", "p_corrected"} <= set(samples.columns)
    assert samples["pvals"].between(0, 1).all()


def test_enrichment_covers_expansion_studies(enriched_search):
    samples = enriched_search.samples
    assert set(samples["series_id"]) == set(enriched_search.expansion_studies["gse_id"])
    assert (samples["geo_accession"] == samples.index).all()


def test_enrichment_columns_are_unique(enriched_search):
    assert enriched_search.samples.columns.is_unique


def test_metadata_only_without_enrichment(metadata_search):
    samples = metadata_search.samples
    assert samples.shape[0] > 0
    assert sum(samples.columns.isin(["extract_protocol_ch1"])) == 1
    assert sum(samples.columns.isin(["p_corrected"])) == 0
    assert set(samples["series_id"]) == set(metadata_search.expansion_studies["gse_id"])


def test_sample_metadata_case(metadata_search):
    # the ARCHS4 reader returns sample metadata as stored (0.2.x upper-cased it)
    assert metadata_search.samples["title"].str.startswith("Sample GSM").all()


def test_samples_tables_have_unnamed_index(enriched_search, metadata_search, tmp_path):
    assert enriched_search.samples.index.name is None
    assert metadata_search.samples.index.name is None
    metadata_search.samples.to_csv(tmp_path / "samples.csv")
    assert (tmp_path / "samples.csv").read_text().startswith(",geo_accession,")


@pytest.fixture
def studies_without_archs4_samples(query_db_archs4, monkeypatch):
    """Make the search select a study that has no samples in the ARCHS4 file."""
    studies = pd.DataFrame({"gse_id": ["GSE_NOT_IN_ARCHS4"]})
    monkeypatch.setattr(query_db_archs4, "transcriptome_search_with_transcriptome_expansion",
                        lambda geneset_query, search=1000, expand=10: (studies, studies))
    return query_db_archs4


def test_enrichment_without_archs4_samples_is_empty(enriched_search, studies_without_archs4_samples, caplog):
    with caplog.at_level(logging.WARNING):
        samples = studies_without_archs4_samples.search(geneset=geneset_query, text_query=None, n_seed=20,
                                                        n_expansion=1, perform_enrichment=True).samples
    assert samples.empty
    assert samples.columns.tolist() == enriched_search.samples.columns.tolist()
    assert "no samples in the ARCHS4 file" in caplog.text


def test_metadata_only_without_archs4_samples_is_empty(studies_without_archs4_samples, caplog):
    with caplog.at_level(logging.WARNING):
        samples = studies_without_archs4_samples.search(geneset=geneset_query, text_query=None, n_seed=20,
                                                        n_expansion=1, perform_enrichment=False).samples
    assert samples.shape == (0, len(SAMPLE_FIELDS))
    assert samples.columns.tolist() == list(SAMPLE_FIELDS)
    assert "no samples in the ARCHS4 file" in caplog.text
