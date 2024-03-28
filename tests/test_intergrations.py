
from biorag.biorag import Query_DB
import os


#query_db = Query_DB("tests/semantic_db.h5ad", "tests/transcriptomic_db.h5ad", "tests/human_gene_v2.2.h5")
text_query = "Studies that investigate severity of COVID-19 in adults"
geneset_query = ["ACE2", "TMPRSS2", "IFITM3", "DPP4", "CTSL", "IL6", "TNF", "IFNG", "CCL2", "STAT1"]

query_db_short = Query_DB("tests/test_semantic_db.h5ad", "tests/test_transcriptome_db.h5ad", "tests/human_gene_v2.2.h5")

#outside
#query_db_short = Query_DB("bioRAG/tests/test_semantic_db.h5ad", "bioRAG/tests/test_transcriptome_db.h5ad", "bioRAG/tests/human_gene_v2.2.h5")


result_COVID = query_db_short.search(geneset = geneset_query, text_query = None, \
                                    n_seed=20, n_expansion=1, perform_enrichment=False)

result_COVID.samples.columns

result_COVID.expansion_studies.shape


def test_search_working_seed_study_dimension():

    result_COVID = query_db_short.search(geneset = geneset_query, text_query = text_query, \
                                   n_seed=20, n_expansion=0, perform_enrichment=True)
    assert result_COVID.seed_studies.shape[0] == 20

def test_search_working_should_be_none_if_n_expansion_none():

    result_COVID = query_db_short.search(geneset = geneset_query, text_query = text_query, \
                                   n_seed=20, n_expansion=0, perform_enrichment=False)
    assert result_COVID.expansion_studies == None

def test_search_working_should_be_none_if_gene_set_none():
    
        result_COVID = query_db_short.search(geneset = None, text_query = text_query, \
                                    n_seed=20, n_expansion=0, perform_enrichment=False)
        assert result_COVID.seed_studies.shape[0] == 20

def test_search_working_should_be_none_if_text_query_none():
    
        result_COVID = query_db_short.search(geneset = None, text_query = text_query, \
                                    n_seed=20, n_expansion=0, perform_enrichment=False)
        assert result_COVID.expansion_studies == None

def test_search_working_should_be_none_if_text_query_none2():
    
        result_COVID = query_db_short.search(geneset = None, text_query = text_query, \
                                    n_seed=20, n_expansion=1, perform_enrichment=False)
        assert result_COVID.expansion_studies.shape[0] == 20

def test_search_working_should_be_none_if_geneset_query_none3():
    
        result_COVID = query_db_short.search(geneset = geneset_query, text_query = None, \
                                    n_seed=20, n_expansion=1, perform_enrichment=False)
        assert result_COVID.expansion_studies.shape[0] == 20

def test_search_working_should_be_none_if_geneset_query_none4():
    
        result_COVID = query_db_short.search(geneset = geneset_query, text_query = None, \
                                    n_seed=20, n_expansion=1, perform_enrichment=False)
        assert result_COVID.expansion_studies.shape[0] == 20 # since n_seed = 20 and expansion would reflect that

def test_search_working_should_be_none_if_geneset_query_none5():
    
        result_COVID = query_db_short.search(geneset = geneset_query, text_query = None, \
                                    n_seed=20, n_expansion=1, perform_enrichment=False)
        assert result_COVID.expansion_studies.shape[0] == 20 # since n_seed = 20 and expansion would reflect that

def test_search_working_should_be_none_if_geneset_query_none6():
    
        result_COVID = query_db_short.search(geneset = geneset_query, text_query = None, \
                                    n_seed=20, n_expansion=1, perform_enrichment=True)
        assert sum(result_COVID.samples.columns.isin(['p_corrected'])) > 0

def test_search_working_should_be_none_if_geneset_query_none7():
    
        result_COVID = query_db_short.search(geneset = geneset_query, text_query = None, \
                                    n_seed=20, n_expansion=1, perform_enrichment=True)
        assert sum(result_COVID.samples.columns.isin(['p_corrected'])) > 0

def test_search_working_should_be_none_if_geneset_query_none8():
    
        result_COVID = query_db_short.search(geneset = geneset_query, text_query = None, \
                                    n_seed=20, n_expansion=1, perform_enrichment=False)
        assert sum(result_COVID.samples.columns.isin(['p_corrected'])) == 0

def test_search_working_should_be_none_if_geneset_query_none9():
    
        result_COVID = query_db_short.search(geneset = geneset_query, text_query = None, \
                                    n_seed=20, n_expansion=1, perform_enrichment=False)
        assert sum(result_COVID.samples.columns.isin(['extract_protocol_ch1'])) == 1