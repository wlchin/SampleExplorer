from pathlib import Path

import anndata as ad
import h5py
import numpy as np

TESTS_DIR = Path(__file__).parent
SEMANTIC_DB = str(TESTS_DIR / "test_semantic_db.h5ad")
TRANSCRIPTOME_DB = str(TESTS_DIR / "test_transcriptome_db.h5ad")

TEST_GENESET = ["IDI1", "SP100", "KLF6", "PLPP1", "NEO1", "TSPAN6"]
TEXT_QUERY = "Trans-chromosomal regulation lincRNA"

ARCHS4_SAMPLE_FIELDS = ["geo_accession", "series_id", "characteristics_ch1",
                        "extract_protocol_ch1", "source_name_ch1", "title"]


def _encode(values):
    return np.array([str(v).encode("utf-8") for v in values])


def build_mini_archs4_h5(path):
    """
    Write a small HDF5 file in the ARCHS4 gene-level layout read by sample_explorer.archs4
    (data/expression genes x samples, meta/genes/symbol, meta/samples/<field>).

    Samples and series come from the transcriptome fixture so that series found by
    a search resolve to samples in this file. Counts are the fixture values for a
    subset of genes, with float16 overflow (inf) capped at the float16 maximum.
    Each sample's most expressed gene is included so that no sample is all zero
    (an all-zero sample has no ranking to test).
    """
    trans = ad.read_h5ad(TRANSCRIPTOME_DB, backed="r")
    all_counts = np.nan_to_num(trans.X[:].astype(np.float32), posinf=65504.0)
    genes = list(trans.var_names[:100]) + list(trans.var_names[990:1000]) + ["IRF6", "IFNG", "CCL2", "STAT1"]
    gene_idx = set(trans.var_names.get_indexer(genes))
    empty = all_counts[:, sorted(gene_idx)].sum(axis=1) == 0
    gene_idx |= set(all_counts[empty].argmax(axis=1))
    gene_idx = sorted(gene_idx)
    genes = list(trans.var_names[gene_idx])
    counts = np.ceil(all_counts[:, gene_idx]).astype(np.uint32)

    obs = trans.obs
    gsm = obs.index.to_list()
    sample_meta = {
        "geo_accession": gsm,
        "series_id": obs["series_id"].to_list(),
        "characteristics_ch1": obs["characteristics_ch1"].to_list(),
        # mixed case on purpose: the old archs4py 0.2.x upper-cased metadata
        "extract_protocol_ch1": [f"Total RNA protocol {i % 7}" for i in range(len(gsm))],
        "source_name_ch1": [f"Source {i % 3}" for i in range(len(gsm))],
        "title": [f"Sample {g}" for g in gsm],
    }

    with h5py.File(path, "w") as f:
        f.create_dataset("data/expression", data=counts.T, compression="gzip", chunks=(len(genes), 100))
        f.create_dataset("meta/genes/symbol", data=_encode(genes))
        for field in ARCHS4_SAMPLE_FIELDS:
            f.create_dataset(f"meta/samples/{field}", data=_encode(sample_meta[field]), compression="gzip")
    trans.file.close()
    return path


def read_archs4_matrix(path, n_samples):
    """Gene symbols, sample IDs and samples x genes counts for the first n_samples, read straight from the HDF5 file."""
    with h5py.File(path, "r") as f:
        genes = [x.decode("utf-8") for x in f["meta/genes/symbol"][:]]
        samples = [x.decode("utf-8") for x in f["meta/samples/geo_accession"][:n_samples]]
        matrix = f["data/expression"][:, :n_samples].T
    return genes, samples, matrix


def read_archs4_counts(path, sample_id):
    """Expression column for one sample, read straight from the HDF5 file."""
    with h5py.File(path, "r") as f:
        gsm = [x.decode("utf-8") for x in f["meta/samples/geo_accession"][:]]
        return f["data/expression"][:, gsm.index(sample_id)]
