from biorag.enrichment import Transcriptome_enrichment
import anndata as ad

x = ad.read_h5ad("tests/test_transcriptome_db.h5ad", backed = "r")

def test_transcriptome_enrichment():
    te = Transcriptome_enrichment(x)
    df = te.list_to_dc_geneset(["IDI1", "SP100","KLF6", "PLPP1", "NEO1", "TSPAN6"])
    df.shape == (6, 2)

gs_dict = {'gsa':["IDI1", "SP100","KLF6", "PLPP1", "NEO1", "TSPAN6"],
           'gsa2':["IDI1", "SP100","KLF6", "PLPP1", "NEO1", "TSPAN6"]}

def test_transcriptome_enrichment_dict():
    te = Transcriptome_enrichment(x)
    df = te.list_to_dc_geneset_dictionary(gs_dict)
    df.shape == (12, 2)