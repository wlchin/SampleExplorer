"""Benchmark sample_explorer.archs4 against archs4py on a file with the ARCHS4 v2.7 layout.

Usage (archs4py is not a project dependency, so add it for the run):

    uv run --with archs4py --with "s3fs>=2024.1" python scripts/benchmark_archs4_reader.py <path.h5>

If <path.h5> does not exist it is built: 67,186 genes x 20,000 samples, uint32 counts with
realistic sparsity, chunks (2000, 1), gzip level 4, plus the six sample metadata fields.
All work is under the __main__ guard because archs4py starts a multiprocessing pool.
"""
import os
import sys
import time

import h5py
import numpy as np

N_GENES = 67_186
N_SAMPLES = 20_000
BLOCK = 250
FIELDS = ["geo_accession", "series_id", "characteristics_ch1",
          "extract_protocol_ch1", "source_name_ch1", "title"]


def build(path):
    rng = np.random.default_rng(0)
    gene_names = [f"GENE{i}" for i in range(N_GENES)]
    gsm = [f"GSM{1000000 + i}" for i in range(N_SAMPLES)]
    meta = {
        "geo_accession": gsm,
        "series_id": [f"GSE{10000 + i // 20}" for i in range(N_SAMPLES)],
        "characteristics_ch1": [f"tissue: t{i % 50}, age: {i % 90}" for i in range(N_SAMPLES)],
        "extract_protocol_ch1": [f"Total RNA protocol {i % 7}" for i in range(N_SAMPLES)],
        "source_name_ch1": [f"Source {i % 30}" for i in range(N_SAMPLES)],
        "title": [f"Sample {g}" for g in gsm],
    }
    with h5py.File(path, "w") as f:
        ds = f.create_dataset("data/expression", shape=(N_GENES, N_SAMPLES), dtype=np.uint32,
                              chunks=(2000, 1), compression="gzip", compression_opts=4)
        base = rng.lognormal(2.0, 2.0, size=N_GENES)
        for start in range(0, N_SAMPLES, BLOCK):
            n = min(BLOCK, N_SAMPLES - start)
            lam = base[:, None] * rng.lognormal(0, 0.3, size=(1, n))
            block = rng.poisson(lam).astype(np.uint32)
            block[rng.random(block.shape) < 0.35] = 0
            ds[:, start:start + n] = block
            if (start // BLOCK) % 10 == 0:
                print(f"  built {start + n}/{N_SAMPLES}", flush=True)
        f.create_dataset("meta/genes/symbol", data=np.array([g.encode() for g in gene_names]))
        for field in FIELDS:
            f.create_dataset(f"meta/samples/{field}", data=np.array([v.encode() for v in meta[field]]))


def main(path):
    import archs4py as a4
    from sample_explorer.archs4 import read_counts, read_sample_metadata

    if not os.path.exists(path):
        print("building", path, flush=True)
        build(path)
    print("file size GB:", round(os.path.getsize(path) / 1e9, 2))

    with h5py.File(path, "r") as f:
        all_ids = [x.decode() for x in f["meta/samples/geo_accession"][:]]
    rng = np.random.default_rng(1)

    print(f"{'n':>6} {'ours (s)':>10} {'archs4py (s)':>13}")
    for n in (100, 1000, 2500):
        ids = list(rng.choice(all_ids, size=n, replace=False))
        t = time.perf_counter()
        ours = read_counts(path, ids)
        t_ours = time.perf_counter() - t
        t = time.perf_counter()
        theirs = a4.data.samples(path, ids, silent=True)
        t_theirs = time.perf_counter() - t
        print(f"{n:>6} {t_ours:>10.2f} {t_theirs:>13.2f}", flush=True)

        assert ours.columns.tolist() == theirs.columns.tolist(), "sample order differs"
        assert ours.index.tolist() == theirs.index.tolist(), "gene symbols differ"
        np.testing.assert_array_equal(ours.to_numpy(), theirs.to_numpy())

        m_ours = read_sample_metadata(path, ids)
        m_theirs = a4.meta.samples(path, ids, silent=True)
        assert m_ours.index.tolist() == m_theirs.index.tolist()
        assert m_ours.columns.tolist() == m_theirs.columns.tolist()
        assert (m_ours.to_numpy() == m_theirs.to_numpy()).all(), "metadata differs"
    print("counts, gene symbols, sample order and metadata identical to archs4py")


if __name__ == "__main__":
    main(sys.argv[1])
