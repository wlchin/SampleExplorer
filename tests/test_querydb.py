

from biorag.biorag import Query_DB
from biorag.biorag import Results

test_geneset =  ["IDI1", "SP100","KLF6", "PLPP1", "NEO1", "TSPAN6"]

#transcriptome_db = ad.read_h5ad("test_transcriptome_db.h5ad", backed = "r")
#te = Transcriptome_enßrichment(transcriptome_db)

# te.memmap_adata.obs.index

# semantic_db = ad.read_h5ad("test_semantic_db.h5ad")
# rag_embedding = Rag_embedding(semantic_db.obs, semantic_db.X)

new_query_db = Query_DB("tests/test_semantic_db.h5ad", "tests/test_transcriptome_db.h5ad")

# res = new_query_db.search(test_geneset, "Trans-chromosomal regulation lincRNA", search = "semantic", expand = "transcriptome", search_n = 50, expand_n = 5)

# res.series
# res.samples

#res.samples.sort_values("query_gene_set", ascending = False)

def test_transcriptome_search():
    res_df, series_of_interest = new_query_db.transcriptome_search(test_geneset, nsamples = 10)
    series_of_interest2 = new_query_db.get_top_samples(res_df, n = 10)
    assert res_df.equals(series_of_interest2)

def test_semantic_search():
    df, res = new_query_db.semantic_search("Trans-chromosomal regulation lincRNA")
    assert df.iloc[0, 2] == "GSE45157"

def test_semantic_search_size_output():
    df, res = new_query_db.semantic_search("Trans-chromosomal regulation lincRNA")
    assert df.shape[0] == 50

def test_transcriptome_expansion_len():
    _, series_of_interest = new_query_db.transcriptome_search(test_geneset, nsamples = 10)
    assert len(series_of_interest) == 10

def test_transcriptome_transcriptome_expansion():
    _, series_of_interest = new_query_db.transcriptome_search(test_geneset, nsamples = 10)
    output = new_query_db.get_transcriptome_series_of_relevance_from_series(series_of_interest, k=5)
    assert output.shape[0] == 50

def test_transcriptome_semantic_expansion():
    _, series_of_interest = new_query_db.transcriptome_search(test_geneset, nsamples = 10)
    output = new_query_db.get_semantic_series_of_relevance_from_series(series_of_interest)
    assert output.shape[0] == 50

def test_semantic_semantic_expansion():
    _, series_of_interest = new_query_db.semantic_search("Trans-chromosomal regulation lincRNA", k = 10)
    output = new_query_db.get_semantic_series_of_relevance_from_series(series_of_interest, k = 5)
    assert output.shape[0] == 50

def test_semantic_transcriptome_expansion():
    _, series_of_interest = new_query_db.semantic_search("Trans-chromosomal regulation lincRNA", k = 10)
    output = new_query_db.get_transcriptome_series_of_relevance_from_series(series_of_interest, k = 5)
    assert output.shape[0] == 50

def test_supporting_transcriptome_transcriptome_expansion():
    _, series_of_interest = new_query_db.transcriptome_search(test_geneset, nsamples = 10)
    output = new_query_db.get_transcriptome_series_of_relevance_from_series(series_of_interest, k=5)
    sup, _ = new_query_db.transcriptome_search_with_transcriptome_expansion(test_geneset, 10, 5)
    assert output.shape[0] == sup.shape[0]

def test_supporting_transcriptome_semantic_expansion():
    _, series_of_interest = new_query_db.transcriptome_search(test_geneset, nsamples = 10)
    output = new_query_db.get_semantic_series_of_relevance_from_series(series_of_interest, 5)
    sup, _ = new_query_db.transcriptome_search_with_semantic_expansion(test_geneset, 10, 5)
    assert output.shape[0] == sup.shape[0]

def test_supporting_semantic_semantic_expansion():
    _, series_of_interest = new_query_db.semantic_search("Trans-chromosomal regulation lincRNA", k = 10)
    output = new_query_db.get_semantic_series_of_relevance_from_series(series_of_interest, k = 5)
    sup, _ = new_query_db.semantic_search_with_semantic_expansion("Trans-chromosomal regulation lincRNA", 10, 5)
    assert output.shape[0] == sup.shape[0]

def test_supporting_semantic_transcriptomic_expansion():
    _, series_of_interest = new_query_db.semantic_search("Trans-chromosomal regulation lincRNA", k = 10)
    output = new_query_db.get_transcriptome_series_of_relevance_from_series(series_of_interest, k = 5)
    sup, _ = new_query_db.semantic_search_with_transcriptome_expansion("Trans-chromosomal regulation lincRNA", 10, 5)
    assert output.shape[0] == sup.shape[0]

def test_search_function():
    """
    Test that the search function returns the same as the semantic search function with expansion
    """
    df = new_query_db.search(test_geneset, "Trans-chromosomal regulation lincRNA", search_n = 5, expand_n = 5)
    df2, _ = new_query_db.semantic_search_with_transcriptome_expansion("Trans-chromosomal regulation lincRNA", 5,5)
    df2 = df2.drop(["similarity_score", "Index"], axis=1).reset_index(drop=True).drop_duplicates() 
    assert df.expansion_studies.equals(df2)

def test_search_function2():
    """
    Test that the search function returns the same as the semantic search function with expansion
    """
    df = new_query_db.search(test_geneset, "Trans-chromosomal regulation lincRNA", search = "semantic", expand = "semantic", search_n = 5, expand_n = 5)
    df2, _ = new_query_db.semantic_search_with_semantic_expansion("Trans-chromosomal regulation lincRNA", 5,5)
    df2 = df2.drop(["similarity_score", "Index"], axis=1).reset_index(drop=True).drop_duplicates() 
    assert df.expansion_studies.equals(df2)

def test_search_function3():
    """
    Test that the search function returns the same as the semantic search function with expansion
    """
    df = new_query_db.search(geneset = test_geneset, text_query="Trans-chromosomal regulation lincRNA", search = "transcriptome", expand = "semantic", search_n = 5, expand_n = 5)
    df2, _ = new_query_db.transcriptome_search_with_semantic_expansion(test_geneset, 5,5)
    df2 = df2.drop(["similarity_score", "Index"], axis=1).reset_index(drop=True).drop_duplicates() 
    assert df.expansion_studies.equals(df2)

def test_search_function4():
    """
    Test that the search function returns the same as the semantic search function with expansion
    """
    df = new_query_db.search(geneset = test_geneset, text_query="Trans-chromosomal regulation lincRNA", search = "transcriptome", expand = "transcriptome", search_n = 5, expand_n = 5)
    df2, _ = new_query_db.transcriptome_search_with_transcriptome_expansion(test_geneset, 5,5)
    df2 = df2.drop(["similarity_score", "Index"], axis=1).reset_index(drop=True).drop_duplicates() 
    assert df.expansion_studies.equals(df2)

def test_results_object():
    """
    Test assignment in Results object
    """

    df = new_query_db.search(geneset = test_geneset, text_query="Trans-chromosomal regulation lincRNA", search = "transcriptome", expand = "transcriptome", search_n = 5, expand_n = 5)
    x = Results(df.seed_studies, df.expansion_studies, df.samples)
    
    assert x.seed_studies.equals(df.seed_studies)