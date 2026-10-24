"""Minimal reader for local ARCHS4 gene-level HDF5 files.

Layout read: ``data/expression`` (genes x samples), ``meta/genes/symbol`` (or
``meta/genes/gene_symbol`` in older files) and ``meta/samples/<field>``.
"""
import json
import os
import shutil
import site
import subprocess
import sys
import tempfile
from collections.abc import Iterable

import h5py
import numpy as np
import pandas as pd

from sample_explorer._archs4_worker import column_runs

SAMPLE_FIELDS = ("geo_accession", "series_id", "characteristics_ch1",
                 "extract_protocol_ch1", "source_name_ch1", "title")

_SAMPLES = "meta/samples/"
# Each sample is a column of gzip chunks, so reads are decompression bound. Below this many
# samples starting worker processes costs more than it saves.
_PARALLEL_MIN_SAMPLES = 500
_MAX_WORKERS = 8


def _decode(values) -> list[str]:
    return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values]


def _read_strings(group_dataset, index=slice(None)) -> list[str]:
    return _decode(group_dataset[index])


def _sample_ids(f) -> np.ndarray:
    return np.array(_read_strings(f[_SAMPLES + "geo_accession"]), dtype=object)


def _positions(f, sample_ids: Iterable[str]) -> tuple[np.ndarray, np.ndarray]:
    """File-order positions of the requested IDs that exist, and all sample IDs."""
    if isinstance(sample_ids, str):
        raise TypeError("sample_ids must be an iterable of sample IDs, not a string.")
    all_ids = _sample_ids(f)
    wanted = set(map(str, sample_ids))
    return np.flatnonzero(np.fromiter((s in wanted for s in all_ids), dtype=bool, count=len(all_ids))), all_ids


def read_series_sample_map(path) -> pd.DataFrame:
    """All samples in file order: index ``geo_accession``, column ``series_id``."""
    with h5py.File(path, "r") as f:
        ids = _read_strings(f[_SAMPLES + "geo_accession"])
        series = _read_strings(f[_SAMPLES + "series_id"])
    return pd.DataFrame({"series_id": series}, index=pd.Index(ids, name="geo_accession"))


def read_sample_metadata(path, sample_ids: Iterable[str], fields=SAMPLE_FIELDS) -> pd.DataFrame:
    """Metadata for the requested samples found in the file, in file order."""
    with h5py.File(path, "r") as f:
        pos, all_ids = _positions(f, sample_ids)
        columns = {}
        for field in fields:
            key = _SAMPLES + field
            if key in f:
                columns[field] = _read_strings(f[key], pos) if len(pos) else []
        return pd.DataFrame(columns, index=pd.Index(all_ids[pos].tolist(), name="geo_accession"))


def samples_in_series(path, series_id: str) -> list[str]:
    """Sample IDs whose series_id equals ``series_id`` exactly, in file order."""
    with h5py.File(path, "r") as f:
        series = np.array(_read_strings(f[_SAMPLES + "series_id"]), dtype=object)
        ids = _read_strings(f[_SAMPLES + "geo_accession"])
    return [ids[i] for i in np.flatnonzero(series == series_id)]


def _read_columns(ds, pos: np.ndarray) -> np.ndarray:
    """The columns at ``pos`` (sorted) as a genes x len(pos) array, one slice per run of adjacent columns."""
    out = np.empty((ds.shape[0], len(pos)), dtype=np.uint32)
    for start, stop, offset in column_runs(pos):
        out[:, offset:offset + stop - start] = ds[:, start:stop]
    return out


def _usable_cpus() -> int:
    """CPUs this process may run on, honouring affinity limits such as SLURM and cgroup cpusets."""
    if hasattr(os, "process_cpu_count"):  # Python 3.13+
        return os.process_cpu_count() or 1
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


def _worker_command(job_path: str) -> list[str]:
    # -P keeps the current directory off the worker's sys.path, so files there cannot shadow its imports.
    return [sys.executable, "-P", "-m", "sample_explorer._archs4_worker", job_path]


def _worker_env() -> dict[str, str]:
    """Environment in which the workers import this copy of sample_explorer, also when the
    caller found it through a sys.path edit rather than an install."""
    env = dict(os.environ)
    root = os.path.realpath(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    site_dirs = {os.path.realpath(p) for p in [*site.getsitepackages(), site.getusersitepackages()]}
    if root not in site_dirs:
        env["PYTHONPATH"] = os.pathsep.join(p for p in (root, env.get("PYTHONPATH")) if p)
    return env


def _tail(path, lines: int = 20) -> str:
    with open(path, "rb") as fh:
        text = fh.read()[-8192:].decode("utf-8", errors="replace")
    return "\n".join(text.splitlines()[-lines:])


def _read_parallel(path, pos: np.ndarray, n_genes: int, workers: int) -> np.ndarray:
    """Read the columns at ``pos`` with ``workers`` subprocesses writing into a temporary file.

    The file is samples x genes, so each worker writes its own contiguous block of it; the
    result is transposed back to genes x samples here.
    """
    shape = (len(pos), n_genes)
    nbytes = n_genes * len(pos) * 4
    tmp = tempfile.mkdtemp(prefix="sample_explorer_archs4_")
    procs = []
    try:
        out_path = os.path.join(tmp, "counts.u32")
        with open(out_path, "wb") as fh:
            fh.truncate(nbytes)
        env = _worker_env()
        bounds = np.linspace(0, len(pos), workers + 1).astype(int)
        for i, (a, b) in enumerate(zip(bounds[:-1], bounds[1:])):
            runs = [(start, stop, offset + int(a)) for start, stop, offset in column_runs(pos[a:b])]
            job_path = os.path.join(tmp, f"job{i}.json")
            with open(job_path, "w", encoding="utf-8") as fh:
                json.dump({"h5": os.path.abspath(path), "out": out_path, "shape": shape, "runs": runs}, fh)
            with open(os.path.join(tmp, f"job{i}.err"), "wb") as err:
                procs.append(subprocess.Popen(_worker_command(job_path), stdin=subprocess.DEVNULL,
                                              stdout=subprocess.DEVNULL, stderr=err, env=env))
        codes = [p.wait() for p in procs]
        failed = [i for i, code in enumerate(codes) if code != 0]
        if failed:
            i, code = failed[0], codes[failed[0]]
            how = f"was killed by signal {-code}" if code < 0 else f"exited with code {code}"
            raise RuntimeError(
                f"Reading {len(pos)} ARCHS4 samples with worker processes failed: {len(failed)} of "
                f"{len(procs)} workers failed, and worker {i} {how}. Likely causes: the system ran out of "
                f"memory; the temporary directory {os.path.dirname(tmp)} ran out of disk space (this read "
                f"needs about {nbytes / 1e6:.1f} MB there; set TMPDIR to use another directory); or the "
                f"ARCHS4 file {path} could not be read.\n"
                f"Worker {i} error output (last lines):\n{_tail(os.path.join(tmp, f'job{i}.err')) or '(none)'}"
            )
        by_sample = np.memmap(out_path, dtype=np.uint32, mode="r", shape=shape)
        counts = np.ascontiguousarray(by_sample.T)
        del by_sample  # unmap before the file is removed
        return counts
    finally:
        for p in procs:
            if p.poll() is None:
                p.kill()
                p.wait()
        shutil.rmtree(tmp, ignore_errors=True)


def read_counts(path, sample_ids: Iterable[str]) -> pd.DataFrame:
    """Counts (genes x found samples, uint32), columns in file order.

    Requests of 500 or more samples are split across worker processes
    (``python -m sample_explorer._archs4_worker``) that write into a temporary file about the
    size of the result in ``tempfile.gettempdir()`` (set TMPDIR to move it).
    """
    with h5py.File(path, "r") as f:
        pos, all_ids = _positions(f, sample_ids)
        symbol_key = "meta/genes/symbol" if "meta/genes/symbol" in f else "meta/genes/gene_symbol"
        genes = _read_strings(f[symbol_key])
        ds = f["data/expression"]
        n_genes = ds.shape[0]
        workers = min(_MAX_WORKERS, _usable_cpus(), len(pos))
        parallel = len(pos) >= _PARALLEL_MIN_SAMPLES and workers > 1 and n_genes > 0
        if not parallel:
            out = _read_columns(ds, pos)
    if parallel:
        out = _read_parallel(path, pos, n_genes, workers)
    return pd.DataFrame(out, index=pd.Index(genes), columns=pd.Index(all_ids[pos].tolist(), name="geo_accession"))
