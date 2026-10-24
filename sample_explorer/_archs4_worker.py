"""Worker process for ``sample_explorer.archs4.read_counts``.

Run as ``python -m sample_explorer._archs4_worker <job.json>``. The job names the ARCHS4 file, a
uint32 file of shape ``shape`` (samples x genes) to write into, and the column runs to copy; each
worker owns a contiguous block of rows of that file. This module imports only numpy and h5py, so
a worker never re-runs the caller's script.
"""
import json
import sys

import h5py
import numpy as np


def column_runs(pos) -> list[tuple[int, int, int]]:
    """Split sorted column positions into runs of adjacent columns: (start, stop, output offset)."""
    runs = []
    start = 0
    while start < len(pos):
        end = start
        while end + 1 < len(pos) and pos[end + 1] == pos[end] + 1:
            end += 1
        runs.append((int(pos[start]), int(pos[end]) + 1, start))
        start = end + 1
    return runs


def main(job_path: str) -> None:
    with open(job_path, encoding="utf-8") as fh:
        job = json.load(fh)
    # No flush: the parent reads these pages back through the same page cache.
    out = np.memmap(job["out"], dtype=np.uint32, mode="r+", shape=tuple(job["shape"]))
    with h5py.File(job["h5"], "r") as f:
        ds = f["data/expression"]
        for start, stop, offset in job["runs"]:
            out[offset:offset + stop - start] = ds[:, start:stop].T
    del out


if __name__ == "__main__":
    main(sys.argv[1])
