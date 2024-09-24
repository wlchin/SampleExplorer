from sample_explorer.sample_explorer import Transcriptome_embedding, Rag_embedding, RNASeqAnalysis
from sample_explorer.utils import MsigDB_store
import pandas as pd
import os
import traceback
import pickle
from typing import List

DATA_DIR = "data/"
LOAD_DIR = "tempfiles/sem_query_rag_df/"
TARGET_DIR = "tempfiles/sem_and_tran_df/"
GENESET_DIR = "gene_sets/"

def load_test_genesets() -> List[str]:
    test_gs = pd.read_csv(f"{GENESET_DIR}test_gene_sets.txt", header=None)[0].to_list()
    return test_gs

def check_folder_existence():
    os.makedirs(TARGET_DIR, exist_ok=True)

def load_rag_data():
    rag_index = pd.read_pickle(f"{DATA_DIR}important_data_for_LLM/rag_index_v2.pkl")
    with open(f"{DATA_DIR}important_data_for_LLM/rag_embedding_matrix_v2.pkl", "rb") as f:
        rag_embedding_matrix = pickle.load(f)
    return Rag_embedding(rag_index, rag_embedding_matrix)

def load_transcriptome_embedding():
    with open("data/important_data_for_LLM/transcription_embedding_matrix_v2.pkl", "rb") as f:
        embedding_matrix = pickle.load(f)

    embedding_index = pd.read_pickle("data/important_data_for_LLM/transcription_index_v2.p")

    exp = Transcriptome_embedding(embedding_index, embedding_matrix)
    return exp

rag = load_rag_data()
msigdb_store = MsigDB_store("data/msigdb.v2023.2.Hs.symbols.gmt", "data/test_set_misgdb.csv")
gene_lists_collection = msigdb_store.list_genesets()
x = RNASeqAnalysis("data/human_gene_v2.2.h5")
exp = load_transcriptome_embedding()
filtered_list = load_test_genesets()
df_series = pd.read_pickle("data/whole_metadata_human.p")
check_folder_existence()

for i in filtered_list:
    try:
        testset = msigdb_store.get_gene_set_by_name(i)        
        basename = i.replace(".csv", "")
        loadfile = f"{LOAD_DIR}{basename}.csv"
        dest_filename = f"{TARGET_DIR}{basename}.csv"
        if os.path.exists(dest_filename):
            print(basename, "exists")
        else:
            rag_df = pd.read_csv(loadfile)
            exp_df = [exp.get_closest_transcriptional_studies(id, k=10) for id in rag_df["gse_id"]]
            studies = pd.concat(exp_df)["gse_id"].drop_duplicates()
            samps = df_series[df_series["series_id"].isin(studies)].index.drop_duplicates()
            if len(samps) < 5000:
                anndata_res = x.perform_enrichment_on_samples(samps, testset["geneset"])
                anndata_res.to_csv(dest_filename)
            else:
                anndata_res = x.perform_enrichment_on_samples_batched(samps, testset["geneset"], batch_size=5000)
                anndata_res.to_csv(dest_filename)
    except:
        basename = filtered_list[i]
        print("failure: ", basename)
        traceback.print_exc()

