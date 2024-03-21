
import anndata as ad
from biorag.biorag import Transcriptome_enrichment


test_geneset =  ["IDI1", "SP100","KLF6", "PLPP1", "NEO1", "TSPAN6"]
x = ad.read_h5ad("tests/test_transcriptome_db.h5ad", backed = "r")
te = Transcriptome_enrichment(x)


def test_dataframe_shape():
    samples_list = te.memmap_adata.obs.head(200).index
    # Replace 'path_to_dataframe' with the actual path to the dataframe produced by your code

    # Load the dataframe
    dataframe = te.run_decouplr_on_memmaped_adata_with_samples(test_geneset, samples_list)

    # Test the shape of the dataframe
    assert dataframe.shape == (200, 3)  # Replace (1000, 2) with the expected shape of the dataframe

def test_observation_anndata_shape():
    obs_df_shape = te.memmap_adata.obs.shape
    # Replace 'path_to_dataframe' with the actual path to the dataframe produced by your code

    assert obs_df_shape == (2000, 2)  # Replace (1000, 2) with the expected shape of the dataframe

def test_gene_list_creation_for_anndata():
    gene_df = te.list_to_dc_geneset(test_geneset)
    # Replace 'path_to_dataframe' with the actual path to the dataframe produced by your code

    assert gene_df.shape == (len(test_geneset), 2)  # Replace (1000, 2) with the expected shape of the dataframe