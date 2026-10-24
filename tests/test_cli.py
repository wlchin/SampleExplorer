import sys

import pandas as pd
import pytest

from sample_explorer import cli
from sample_explorer.sample_explorer import Results
from tests.common import SEMANTIC_DB, TEST_GENESET, TEXT_QUERY, TRANSCRIPTOME_DB


class FakeQueryDB:
    """Records how the CLI builds and calls Query_DB, without loading any stores."""
    instances = []

    def __init__(self, semantic_db, transcriptome_db, archs4_file=None):
        self.args = (semantic_db, transcriptome_db, archs4_file)
        self.search_kwargs = None
        FakeQueryDB.instances.append(self)

    def search(self, **kwargs):
        self.search_kwargs = kwargs
        df = pd.DataFrame({"gse_id": ["GSE1"]})
        return Results(df, df, df)


@pytest.fixture
def fake_query_db(monkeypatch):
    FakeQueryDB.instances = []
    monkeypatch.setattr(cli, "Query_DB", FakeQueryDB)
    return FakeQueryDB


def run_cli(monkeypatch, tmp_path, *args):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["biorag", *args])
    cli.main()
    (output_folder,) = (tmp_path / "results").iterdir()
    return output_folder


def test_usage_exits_cleanly(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["biorag", "--usage"])
    with pytest.raises(SystemExit) as exit_info:
        cli.main()
    assert exit_info.value.code == 0
    assert "Usage on command line" in capsys.readouterr().out


def test_comma_separated_gene_list_and_text(monkeypatch, tmp_path, fake_query_db):
    run_cli(monkeypatch, tmp_path, "--gene_list", "IRF1,IRF2,IRF3", "--text_query", "cardiac myocytes",
            "--semantic_db", "sem.h5ad", "--transcriptome_db", "trans.h5ad", "--archs4_file", "archs4.h5")
    (db,) = fake_query_db.instances
    assert db.args == ("sem.h5ad", "trans.h5ad", "archs4.h5")
    assert db.search_kwargs == {"geneset": ["IRF1", "IRF2", "IRF3"], "text_query": "cardiac myocytes",
                                "search": "semantic", "expand": "transcriptome", "perform_enrichment": False}


def test_gene_list_and_text_query_files(monkeypatch, tmp_path, fake_query_db):
    (tmp_path / "genes.txt").write_text("IRF1\nIRF2 \n")
    (tmp_path / "query.txt").write_text("studies with cardiac myocytes\n")
    run_cli(monkeypatch, tmp_path, "--gene_list", "genes.txt", "--text_query", "query.txt",
            "--search", "transcriptome", "--expand", "semantic", "--enrichment")
    (db,) = fake_query_db.instances
    assert db.search_kwargs["geneset"] == ["IRF1", "IRF2"]
    assert db.search_kwargs["text_query"] == "studies with cardiac myocytes\n"
    assert db.search_kwargs["search"] == "transcriptome"
    assert db.search_kwargs["expand"] == "semantic"
    assert db.search_kwargs["perform_enrichment"] is True


def test_writes_result_tables(monkeypatch, tmp_path, fake_query_db):
    output_folder = run_cli(monkeypatch, tmp_path, "--gene_list", "IRF1,IRF2", "--text_query", "x")
    assert sorted(p.name for p in output_folder.iterdir()) == [
        "relevant_expansion_studies.csv", "relevant_samples.csv", "relevant_seed_studies.csv"]


def test_end_to_end_with_archs4_file(monkeypatch, tmp_path, archs4_h5):
    output_folder = run_cli(monkeypatch, tmp_path, "--gene_list", ",".join(TEST_GENESET), "--text_query", TEXT_QUERY,
                            "--semantic_db", SEMANTIC_DB, "--transcriptome_db", TRANSCRIPTOME_DB,
                            "--archs4_file", archs4_h5)
    seeds = pd.read_csv(output_folder / "relevant_seed_studies.csv")
    expansion = pd.read_csv(output_folder / "relevant_expansion_studies.csv")
    samples = pd.read_csv(output_folder / "relevant_samples.csv")
    # defaults: 50 semantic seed studies, each expanded to 10 transcriptome neighbours
    assert len(seeds) == 50
    assert len(expansion) == 500
    assert set(samples["series_id"]) == set(expansion["gse_id"])


@pytest.mark.xfail(raises=AttributeError, reason="without --archs4_file, Results.samples is None and main() calls None.to_csv")
def test_end_to_end_without_archs4_file(monkeypatch, tmp_path):
    output_folder = run_cli(monkeypatch, tmp_path, "--gene_list", ",".join(TEST_GENESET), "--text_query", TEXT_QUERY,
                            "--semantic_db", SEMANTIC_DB, "--transcriptome_db", TRANSCRIPTOME_DB)
    assert (output_folder / "relevant_seed_studies.csv").exists()
