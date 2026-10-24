"""Minimal reader for local ARCHS4 gene-level HDF5 files.

Layout read: ``data/expression`` (genes x samples), ``meta/genes/symbol`` (or
``meta/genes/gene_symbol`` in older files) and ``meta/samples/<field>``.
"""
import multiprocessing as mp
import os
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import shared_memory

import h5py
import numpy as np
import pandas as pd

SAMPLE_FIELDS = ("geo_accession", "series_id", "characteristics_ch1",
                 "extract_protocol_ch1", "source_name_ch1", "title")

_SAMPLES = "meta/samples/"
# Each sample is a column of gzip chunks, so reads are decompression bound. Below this many
# samples a process pool costs more than it saves.
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


def _read_columns(ds, pos: np.ndarray, out: np.ndarray, offset: int = 0) -> None:
    """Fill out[:, offset:offset+len(pos)] with the columns at ``pos`` (sorted), one slice per run."""
    start = 0
    while start < len(pos):
        end = start
        while end + 1 < len(pos) and pos[end + 1] == pos[end] + 1:
            end += 1
        out[:, offset + start:offset + end + 1] = ds[:, pos[start]:pos[end] + 1]
        start = end + 1


def _worker(args) -> None:
    path, pos, shm_name, shape, offset = args
    shm = shared_memory.SharedMemory(name=shm_name)
    try:
        out = np.ndarray(shape, dtype=np.uint32, buffer=shm.buf)
        with h5py.File(path, "r") as f:
            _read_columns(f["data/expression"], pos, out, offset)
    finally:
        shm.close()


def _read_parallel(path, pos: np.ndarray, n_genes: int, workers: int) -> np.ndarray:
    shape = (n_genes, len(pos))
    shm = shared_memory.SharedMemory(create=True, size=max(1, int(np.prod(shape)) * 4))
    try:
        bounds = np.linspace(0, len(pos), workers + 1).astype(int)
        jobs = [(path, pos[a:b], shm.name, shape, int(a)) for a, b in zip(bounds[:-1], bounds[1:]) if b > a]
        # spawn: safe with threaded libraries loaded; callers need an `if __name__ == "__main__":` guard
        with ProcessPoolExecutor(len(jobs), mp_context=mp.get_context("spawn")) as pool:
            list(pool.map(_worker, jobs))
        return np.ndarray(shape, dtype=np.uint32, buffer=shm.buf).copy()
    finally:
        shm.close()
        shm.unlink()


def read_counts(path, sample_ids: Iterable[str]) -> pd.DataFrame:
    """Counts (genes x found samples, uint32), columns in file order.

    Large requests are read by a spawn process pool, so scripts calling this with 500 or more
    samples need an ``if __name__ == "__main__":`` guard.
    """
    with h5py.File(path, "r") as f:
        pos, all_ids = _positions(f, sample_ids)
        symbol_key = "meta/genes/symbol" if "meta/genes/symbol" in f else "meta/genes/gene_symbol"
        genes = _read_strings(f[symbol_key])
        ds = f["data/expression"]
        workers = min(_MAX_WORKERS, os.cpu_count() or 1)
        if len(pos) >= _PARALLEL_MIN_SAMPLES and workers > 1:
            out = None
            n_genes = ds.shape[0]
        else:
            out = np.empty((ds.shape[0], len(pos)), dtype=np.uint32)
            _read_columns(ds, pos, out)
    if out is None:
        out = _read_parallel(path, pos, n_genes, workers)
    return pd.DataFrame(out, index=pd.Index(genes), columns=pd.Index(all_ids[pos].tolist(), name="geo_accession"))
