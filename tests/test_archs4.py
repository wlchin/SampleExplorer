import logging

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest

from sample_explorer.archs4 import (SAMPLE_FIELDS, read_counts, read_sample_metadata,
                                    read_series_sample_map, samples_in_series)
from sample_explorer.rnaseq_analysis import RNASeqAnalysis
from tests.common import TRANSCRIPTOME_DB, _encode, read_archs4_counts


def test_series_sample_map(archs4_h5):
    obs = ad.read_h5ad(TRANSCRIPTOME_DB, backed="r").obs
    result = read_series_sample_map(archs4_h5)
    assert len(result) == 2000
    assert result.index.tolist() == obs.index.tolist()
    assert result.columns.tolist() == ["series_id"]
    assert result["series_id"].tolist() == obs["series_id"].tolist()


def test_sample_metadata_order_values_and_unknown_ids(archs4_h5):
    meta = read_sample_metadata(archs4_h5, pd.Index(["GSM1002540", "GSM1000981", "GSM_UNKNOWN"]))
    assert meta.index.tolist() == ["GSM1000981", "GSM1002540"]
    assert meta.columns.tolist() == list(SAMPLE_FIELDS)
    assert meta.loc["GSM1000981", "title"] == "Sample GSM1000981"


def test_samples_in_series(archs4_h5):
    assert samples_in_series(archs4_h5, "GSE29282") == ["GSM1000981", "GSM1000984"]
    assert samples_in_series(archs4_h5, "GSE_NONE") == []


def test_read_counts_matches_file(archs4_h5):
    ids = np.array(["GSM1002540", "GSM1000984", "GSM1000981"])
    counts = read_counts(archs4_h5, ids)
    assert counts.columns.tolist() == ["GSM1000981", "GSM1000984", "GSM1002540"]
    assert all(counts.dtypes == np.uint32)
    for sample_id in counts.columns:
        np.testing.assert_array_equal(counts[sample_id].to_numpy(), read_archs4_counts(archs4_h5, sample_id))


def test_older_layout_and_missing_fields(tmp_path):
    path = tmp_path / "old.h5"
    gsm = ["GSM1", "GSM2", "GSM3"]
    with h5py.File(path, "w") as f:
        f.create_dataset("data/expression", data=np.array([[5, 1, 0], [3, 2, 9], [1, 1, 4]], dtype=np.uint32))
        f.create_dataset("meta/genes/gene_symbol", data=_encode(["A", "A", "B"]))
        for field in ("geo_accession", "series_id", "characteristics_ch1", "source_name_ch1", "title"):
            values = gsm if field == "geo_accession" else [f"{field} {g}" for g in gsm]
            f.create_dataset(f"meta/samples/{field}", data=_encode(values))
    counts = read_counts(str(path), gsm)
    assert counts.index.tolist() == ["A", "A", "B"]
    meta = read_sample_metadata(str(path), gsm)
    assert "extract_protocol_ch1" not in meta.columns
    assert "title" in meta.columns
    result = RNASeqAnalysis(str(path)).perform_enrichment_on_samples(gsm, ["A", "B"])
    assert len(result) == 3


def test_create_anndata_warns_and_raises(archs4_h5, caplog):
    rnaseq = RNASeqAnalysis(archs4_h5)
    with caplog.at_level(logging.WARNING):
        adata = rnaseq.create_anndata_from_samples(["GSM1000981", "GSM_UNKNOWN"])
    assert adata.n_obs == 1
    assert any("1" in r.message and "not found" in r.message for r in caplog.records)
    with pytest.raises(ValueError):
        rnaseq.create_anndata_from_samples(["GSM_UNKNOWN"])


def test_parallel_read_matches_serial(archs4_h5, monkeypatch):
    ids = ["GSM1002540", "GSM1000984", "GSM1000981", "GSM1000985"]
    serial = read_counts(archs4_h5, ids)
    monkeypatch.setattr("sample_explorer.archs4._PARALLEL_MIN_SAMPLES", 1)
    parallel = read_counts(archs4_h5, ids)
    pd.testing.assert_frame_equal(serial, parallel)
