import pytest
import pandas as pd
from biorag import Rag_embedding, Transcriptome_enrichment

@pytest.fixture
def rag_embedding_instance():
    # Set up the Rag_embedding instance with appropriate paths
    return Rag_embedding("workflow/data/important_data_for_LLM/rag_index.p", 
                    "workflow/data/important_data_for_LLM/rag_embedding_matrix.pkl")

@pytest.fixture
def transcriptome_enrichment_instance():
    # Set up the Rag_embedding instance with appropriate paths
    return Transcriptome_enrichment("workflow/data/study_counts.dat", 
                         "workflow/data/important_data_for_LLM/genes.p", 
                         "workflow/data/important_data_for_LLM/transcription_index.p", "workflow/data/smaller.h5ad")


# def test_get_gse_from_rag_index(rag_embedding_instance):
#     # Test if the method returns the correct GSE ID
#     gse_id = rag_embedding_instance.get_gse_from_rag_index(0)
#     assert gse_id == expected_gse_id

# def test_get_gse_text_from_rag_index(rag_embedding_instance):
#     # Test if the method returns the correct GSE ID text
#     gse_id_text = rag_embedding_instance.get_gse_text_from_rag_index(0)
#     assert gse_id_text == expected_gse_id_text

# def test_get_text_linked_to_gse(rag_embedding_instance):
#     # Test if the method returns the correct text linked to GSE ID
#     gse_id_query = 'your_gse_id_here'
#     text_linked = rag_embedding_instance.get_text_linked_to_gse(gse_id_query)
#     assert text_linked == expected_text_linked

def test_get_closest_semantic_studies(rag_embedding_instance):
    # Test if the method returns the expected DataFrame structure
    query_gse_id = 'GSE90133'
    result_df = rag_embedding_instance.get_closest_semantic_studies(query_gse_id, k=5)
    assert isinstance(result_df, pd.DataFrame)
    # Add more assertions based on the expected structure of the result DataFrame
    
test_dict = {"genelistA": ["TRNE", "CYTB", "TRNT", "TRNP", "DDX11L1"], 
             "genelistB": ["ND6", "CYTB", "TRNT", "TRNP", "DDX11L1"]}

test_list = ["ND6", "TRNE", "CYTB", "TRNT", "TRNP", "DDX11L1"]

series_list_test = ["GSE90133"]
samples_list_test = ["GSM1958505", "GSM1958508"]

def test_mutligeneset_construction(transcriptome_enrichment_instance):
    # Test if the method returns the expected DataFrame structure
    result_df = transcriptome_enrichment_instance.list_to_dc_geneset_dictionary(test_dict)
    assert isinstance(result_df, pd.DataFrame)
    # Add more assertions based on the expected structure of the result DataFrame
    
def test_mutligeneset_construction_geneset_dim(transcriptome_enrichment_instance):
    # Test if the method returns the expected DataFrame structure
    result_df = transcriptome_enrichment_instance.list_to_dc_geneset_dictionary(test_dict)
    assert result_df.shape == (10,2)
    # Add more assertions based on the expected structure of the result DataFrame
    
def test_enrichment_list_format(transcriptome_enrichment_instance):
    # Test if the method returns the expected DataFrame structure
    result_df = transcriptome_enrichment_instance.run_decouplr_on_memmaped_adata(test_list, series_list_test)
    assert result_df.shape == (12,3)
    # Add more assertions based on the expected structure of the result DataFrame
    
def test_enrichment_dict_formats(transcriptome_enrichment_instance):
    # Test if the method returns the expected DataFrame structure
    res = transcriptome_enrichment_instance.run_decouplr_on_memmaped_adata(test_dict, series_list_test)
    assert res[0].shape == (12, 2)
    # Add more assertions based on the expected structure of the result DataFrame
    
def test_enrichment_samples_genelist_and_samples(transcriptome_enrichment_instance):
    # Test if the method returns the expected DataFrame structure
    res = transcriptome_enrichment_instance.run_decouplr_on_memmaped_adata_with_samples(test_list, samples_list_test)
    assert res.shape == (2, 3)
    # Add more assertions based on the expected structure of the result DataFrame
    
def test_enrichment_samples_genedict_and_samples(transcriptome_enrichment_instance):
    # Test if the method returns the expected DataFrame structure
    res = transcriptome_enrichment_instance.run_decouplr_on_memmaped_adata_with_samples(test_dict, samples_list_test)
    assert res[0].shape == (2, 2)
    # Add more assertions based on the expected structure of the result DataFrame
    
