


from biorag.biorag import Rag_embedding
import anndata as ad

x = ad.read_h5ad("tests/test_semantic_db.h5ad")

rag_embedding = Rag_embedding(x.obs, x.X)

def test_get_gse_from_rag_index():
    ans = rag_embedding.get_gse_from_rag_index(3)
    assert ans == 'GSE44615,GSE44616'

def test_rag_returned():
    gse_id = rag_embedding.get_gse_from_rag_index(3)
    text1 = rag_embedding.get_text_linked_to_gse(gse_id)
    text2 = rag_embedding.get_gse_text_from_rag_index(3)
    assert text1 == text2  # Replace (1000, 2) with the expected shape of the dataframe

def test_rag_sem_df_size():
    df = rag_embedding.get_closest_semantic_studies('GSE44615,GSE44616', 7)
    assert df.shape == (7, 4)

def test_query_rag():
    df = rag_embedding.query_rag("Trans-chromosomal regulation lincRNA")
    assert df.iloc[0, 2] == "GSE45157"