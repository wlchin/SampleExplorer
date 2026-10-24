import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tests.common import TEST_GENESET


def test_transcriptome_search_propagates_errors(query_db, monkeypatch):
    def boom(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(
        query_db.transcriptome_enrichment,
        "run_decouplr_on_memmaped_adata_with_samples",
        boom,
    )
    with pytest.raises(RuntimeError, match="boom"):
        query_db.transcriptome_search(TEST_GENESET, nsamples=5)


def test_semantic_expansion_skips_unknown_series(query_db, caplog):
    out = query_db.get_semantic_series_of_relevance_from_series(["GSE45157", "GSE_NOT_THERE"], k=3)
    assert len(out) == 3 and "GSE_NOT_THERE" in caplog.text


def test_transcriptome_expansion_skips_unknown_series(query_db, caplog):
    out = query_db.get_transcriptome_series_of_relevance_from_series(["GSE43801", "GSE_NOT_THERE"], k=3)
    assert len(out) == 3 and "GSE_NOT_THERE" in caplog.text


def test_expansion_of_only_unknown_series_keeps_columns(query_db):
    sem = query_db.get_semantic_series_of_relevance_from_series(["GSE_NOT_THERE"], k=3)
    assert sem.empty
    assert list(sem.columns) == ["Index", "similarity_score", "gse_id", "gse_id_text"]
    tr = query_db.get_transcriptome_series_of_relevance_from_series(["GSE_NOT_THERE"], k=3)
    assert tr.empty
    assert list(tr.columns) == ["Index", "similarity_score", "gse_id", "gse_id_text"]


def test_semantic_expansion_clamps_k(query_db):
    assert len(query_db.rag_embedding.get_closest_semantic_studies("GSE45157", k=500)) == 98


def test_text_lookup_unknown_series_returns_none(query_db):
    assert query_db.rag_embedding.get_text_linked_to_gse("GSE_NOT_THERE") is None


def test_top_samples_ties_keep_store_order(query_db):
    index = query_db.transcriptome_embedding.embeddings_index.index[:5]
    ties = pd.DataFrame({"pvals": 1.0, "estimate": 0.0}, index=index)
    assert query_db.get_top_samples(ties, n=2).index.tolist() == index[:2].tolist()


def _seed_studies(query_db, geneset):
    return query_db.search(geneset=geneset, text_query=None, search="transcriptome", n_seed=5, n_expansion=0).seed_studies


def test_search_accepts_any_gene_iterable(query_db):
    expected = _seed_studies(query_db, list(TEST_GENESET))
    for genes in (tuple(TEST_GENESET), np.array(TEST_GENESET), pd.Series(TEST_GENESET)):
        pd.testing.assert_frame_equal(_seed_studies(query_db, genes), expected)


def test_search_rejects_string_geneset(query_db):
    with pytest.raises(TypeError):
        query_db.search(geneset="IRF1", text_query=None, n_seed=5, n_expansion=0)


def test_search_deduplicates_genes(query_db):
    expected = _seed_studies(query_db, list(TEST_GENESET))
    pd.testing.assert_frame_equal(_seed_studies(query_db, TEST_GENESET + TEST_GENESET[:2]), expected)


def test_import_leaves_warnings_and_logging_alone():
    code = (
        "import logging, warnings\n"
        "root = logging.getLogger()\n"
        "before = (list(root.handlers), root.level)\n"
        "import sample_explorer.sample_explorer\n"
        "assert ('ignore', None, Warning, None, 0) not in warnings.filters, 'blanket ignore filter'\n"
        "assert (list(root.handlers), root.level) == before, 'root logger changed'\n"
    )
    repo_root = Path(__file__).resolve().parent.parent
    result = subprocess.run([sys.executable, "-c", code], cwd=repo_root, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_query_db_adds_no_handlers(query_db):
    assert logging.getLogger("sample_explorer.sample_explorer").handlers == []
    assert logging.getLogger("sample_explorer.rag_embedding").handlers == []
