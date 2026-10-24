import pandas as pd

from sample_explorer.utils import ARCHS4_API_query, MsigDB_store


def write_msigdb_files(tmp_path):
    gmt = tmp_path / "sets.gmt"
    gmt.write_text("SET_A\thttps://example.org/SET_A\tIRF1\tIRF2\tIRF3\n"
                   "SET_B\thttps://example.org/SET_B\tSTAT1\tSTAT2\n")
    meta = tmp_path / "sets.csv"
    pd.DataFrame({
        "standard_name": ["SET_A", "SET_B"],
        "systematic_name": ["M1", "M2"],
        "collection": ["H", "H"],
        "source": ["x", "y"],
        "n_genes": [3, 2],
        "description_full": ["Long title A", "Long title B"],
        "description_brief": ["Short A", "Short B"],
    }).to_csv(meta, index=False)
    return gmt, meta


def test_msigdb_store(tmp_path):
    gmt, meta = write_msigdb_files(tmp_path)
    store = MsigDB_store(gmt, meta)
    assert store.list_genesets().tolist() == ["SET_A", "SET_B"]
    assert store.get_gene_set_by_name("SET_B") == {
        "geneset": ["STAT1", "STAT2"], "short_title": "Short B", "long_title": "Long title B"}


def test_archs4_api_query_metadata(tmp_path):
    pickle_path = tmp_path / "meta.p"
    pd.DataFrame({"series_id": ["GSE1", "GSE1", "GSE2"], "title": ["a", "b", "c"]},
                 index=["GSM1", "GSM2", "GSM3"]).to_pickle(pickle_path)
    query = ARCHS4_API_query(pickle_path)
    assert query.geo_metadata.columns.tolist() == ["series_id", "samples"]
    assert query.geo_metadata["samples"].tolist() == ["GSM1", "GSM2", "GSM3"]
    assert query.create_data_dict("up", ["IRF1"])["upgenes"] == ["IRF1"]
    assert query.create_data_dict("down", ["IRF1"])["downgenes"] == ["IRF1"]
