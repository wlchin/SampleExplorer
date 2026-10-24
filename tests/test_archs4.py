import logging
import os
import subprocess
import sys
import tempfile
import textwrap
import warnings

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest

from sample_explorer import archs4
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


def _old_layout_file(tmp_path):
    """Older layout (meta/genes/gene_symbol), a duplicated symbol and no extract_protocol_ch1."""
    path = tmp_path / "old.h5"
    gsm = ["GSM1", "GSM2", "GSM3"]
    with h5py.File(path, "w") as f:
        f.create_dataset("data/expression", data=np.array([[5, 1, 0], [3, 2, 9], [1, 1, 4]], dtype=np.uint32))
        f.create_dataset("meta/genes/gene_symbol", data=_encode(["A", "A", "B"]))
        for field in ("geo_accession", "series_id", "characteristics_ch1", "source_name_ch1", "title"):
            values = gsm if field == "geo_accession" else [f"{field} {g}" for g in gsm]
            f.create_dataset(f"meta/samples/{field}", data=_encode(values))
    return str(path), gsm


def test_older_layout_and_missing_fields(tmp_path):
    path, gsm = _old_layout_file(tmp_path)
    counts = read_counts(path, gsm)
    assert counts.index.tolist() == ["A", "A", "B"]
    meta = read_sample_metadata(path, gsm)
    assert "extract_protocol_ch1" not in meta.columns
    assert "title" in meta.columns
    result = RNASeqAnalysis(path).perform_enrichment_on_samples(gsm, ["A", "B"])
    assert len(result) == 3


def test_duplicate_gene_symbols_made_unique_without_warning(tmp_path):
    path, gsm = _old_layout_file(tmp_path)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        adata = RNASeqAnalysis(path).create_anndata_from_samples(gsm)
    assert not [w for w in caught if "not unique" in str(w.message)]
    # the names AnnData.var_names_make_unique gives
    assert adata.var_names.tolist() == ["A", "A-1", "B"]


def test_create_anndata_warns_and_raises(archs4_h5, caplog):
    rnaseq = RNASeqAnalysis(archs4_h5)
    with caplog.at_level(logging.WARNING):
        adata = rnaseq.create_anndata_from_samples(["GSM1000981", "GSM_UNKNOWN"])
    assert adata.n_obs == 1
    assert any("1" in r.message and "not found" in r.message for r in caplog.records)
    with pytest.raises(ValueError):
        rnaseq.create_anndata_from_samples(["GSM_UNKNOWN"])


def _all_ids(path):
    with h5py.File(path, "r") as f:
        return [x.decode("utf-8") for x in f["meta/samples/geo_accession"][:]]


def _force_workers(monkeypatch, cpus=4):
    """Use worker processes for any request size, also on a machine with one CPU."""
    monkeypatch.setattr(archs4, "_PARALLEL_MIN_SAMPLES", 1)
    monkeypatch.setattr(archs4, "_usable_cpus", lambda: cpus)


def test_parallel_read_matches_serial(archs4_h5, monkeypatch, tmp_path):
    all_ids = _all_ids(archs4_h5)
    # runs of adjacent columns and scattered columns, requested out of file order
    ids = all_ids[1500:1800] + all_ids[3:1400:7]
    ids.reverse()
    serial = read_counts(archs4_h5, ids)
    _force_workers(monkeypatch)
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    calls = []
    real_read_parallel = archs4._read_parallel
    monkeypatch.setattr(archs4, "_read_parallel", lambda *args: calls.append(args[3]) or real_read_parallel(*args))
    parallel = read_counts(archs4_h5, ids)
    assert calls == [4]
    pd.testing.assert_frame_equal(serial, parallel)
    assert list(tmp_path.iterdir()) == []


def test_failed_worker_raises_with_its_error_output(archs4_h5, monkeypatch, tmp_path):
    _force_workers(monkeypatch, cpus=2)
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    real_command = archs4._worker_command
    fail = "import sys; sys.stderr.write('OSError: Unable to open file (worker test)\\n'); sys.exit(3)"

    def command(job_path):
        return [sys.executable, "-c", fail] if job_path.endswith("job0.json") else real_command(job_path)

    monkeypatch.setattr(archs4, "_worker_command", command)
    with pytest.raises(RuntimeError) as excinfo:
        read_counts(archs4_h5, _all_ids(archs4_h5)[:10])
    message = str(excinfo.value)
    assert "OSError: Unable to open file (worker test)" in message
    assert "1 of 2 workers failed" in message and "exited with code 3" in message
    assert "disk space" in message and "memory" in message
    assert "__main__" not in message
    assert list(tmp_path.iterdir()) == []


def test_unguarded_script_is_not_rerun_by_workers(archs4_h5, tmp_path):
    ids = _all_ids(archs4_h5)[:50]
    counter = tmp_path / "runs.txt"
    script = tmp_path / "unguarded.py"
    script.write_text(textwrap.dedent(f"""\
        # deliberately no `if __name__ == "__main__":` guard
        from sample_explorer import archs4

        with open({str(counter)!r}, "a") as fh:
            fh.write("ran\\n")
        archs4._PARALLEL_MIN_SAMPLES = 1
        archs4._usable_cpus = lambda: 4
        calls = []
        real_read_parallel = archs4._read_parallel
        archs4._read_parallel = lambda *args: calls.append(1) or real_read_parallel(*args)
        counts = archs4.read_counts({archs4_h5!r}, {ids!r})
        assert calls and counts.shape[1] == {len(ids)}, (calls, counts.shape)
        """))
    result = subprocess.run([sys.executable, str(script)], cwd=tmp_path, capture_output=True, text=True,
                            timeout=300)
    assert result.returncode == 0, result.stderr
    assert counter.read_text().splitlines() == ["ran"]


def test_string_sample_id_is_rejected(archs4_h5):
    with pytest.raises(TypeError):
        read_counts(archs4_h5, "GSM1000981")
    with pytest.raises(TypeError):
        read_sample_metadata(archs4_h5, "GSM1000981")
    with pytest.raises(TypeError):
        RNASeqAnalysis(archs4_h5).create_anndata_from_samples("GSM1000981")


def test_any_sample_id_iterable_is_accepted(archs4_h5):
    ids = ["GSM1002540", "GSM1000981"]
    expected = read_counts(archs4_h5, ids)
    for other in (tuple(ids), pd.Index(ids), (i for i in ids)):
        pd.testing.assert_frame_equal(read_counts(archs4_h5, other), expected)
    adata = RNASeqAnalysis(archs4_h5).create_anndata_from_samples(i for i in ids)
    assert adata.obs_names.tolist() == ["GSM1000981", "GSM1002540"]


def test_usable_cpus_honours_affinity(monkeypatch):
    monkeypatch.setattr(os, "cpu_count", lambda: 64)
    monkeypatch.setattr(os, "process_cpu_count", lambda: 3, raising=False)
    assert archs4._usable_cpus() == 3
    monkeypatch.delattr(os, "process_cpu_count")
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0, 1}, raising=False)
    assert archs4._usable_cpus() == 2
    monkeypatch.delattr(os, "sched_getaffinity")
    assert archs4._usable_cpus() == 64
