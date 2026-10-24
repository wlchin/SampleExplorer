from sample_explorer.enrichment import Transcriptome_enrichment
import anndata as ad
from tests.common import TRANSCRIPTOME_DB

x = ad.read_h5ad(TRANSCRIPTOME_DB, backed = "r")

def test_transcriptome_enrichment():
    te = Transcriptome_enrichment(x)
    df = te.list_to_dc_geneset(["IDI1", "SP100","KLF6", "PLPP1", "NEO1", "TSPAN6"])
    assert df.shape == (6, 2)
    assert set(df["geneset"]) == {"query_gene_set"}

gs_dict = {'gsa':["IDI1", "SP100","KLF6", "PLPP1", "NEO1", "TSPAN6"],
           'gsa2':["IDI1", "SP100","KLF6", "PLPP1", "NEO1", "TSPAN6"]}

def test_transcriptome_enrichment_dict():
    te = Transcriptome_enrichment(x)
    df = te.list_to_dc_geneset_dictionary(gs_dict)
    assert df.shape == (12, 2)
    assert df.groupby("geneset").size().to_dict() == {"gsa": 6, "gsa2": 6}
