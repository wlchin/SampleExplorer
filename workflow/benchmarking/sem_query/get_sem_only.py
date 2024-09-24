from biorag.biorag import Rag_embedding, RNASeqAnalysis
from biorag.utils import MsigDB_store
import pandas as pd
import traceback
import pandas as pd
import os
from typing import List
import pickle

# Constants
DATA_DIR = "data"
LOAD_DIR = "tempfiles/sem_query_rag_df/"
TARGET_DIR = "tempfiles/sem_only_df/"
GENESET_DIR = "gene_sets/"

def load_test_genesets() -> List[str]:
    test_gs = pd.read_csv(f"{GENESET_DIR}test_gene_sets.txt", header=None)[0].to_list()
    return test_gs

def check_folder_existence():
    os.makedirs(TARGET_DIR, exist_ok=True)

def load_rag_data():
    rag_index = pd.read_pickle(f"{DATA_DIR}/important_data_for_LLM/rag_index_v2.pkl")
    with open(f"{DATA_DIR}/important_data_for_LLM/rag_embedding_matrix_v2.pkl", "rb") as f:
        rag_embedding_matrix = pickle.load(f)
    return Rag_embedding(rag_index, rag_embedding_matrix)

rag = load_rag_data()

msigdb_store = MsigDB_store("data/msigdb.v2023.2.Hs.symbols.gmt", "data/test_set_misgdb.csv")

x = RNASeqAnalysis("data/human_gene_v2.2.h5")

df_series = pd.read_pickle("data/whole_metadata_human.p")

filtered_list = load_test_genesets()

for i in filtered_list :
    try:
        check_folder_existence()
        testset = msigdb_store.get_gene_set_by_name(i)        
        basename = i.replace(".csv", "")
        loadfile = f"{LOAD_DIR}{basename}.csv"
        dest_filename = f"{TARGET_DIR}{basename}.csv" 
        rag_df = pd.read_csv(loadfile)
        samps = df_series[df_series["series_id"].isin(rag_df["gse_id"])].index.drop_duplicates()
        anndata_res = x.perform_enrichment_on_samples(samps, testset["geneset"])
        anndata_res.to_csv(dest_filename)
    except:
        basename = filtered_list[i]
        print("failure: ", basename)
