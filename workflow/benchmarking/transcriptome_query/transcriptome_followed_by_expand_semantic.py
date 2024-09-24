from sample_explorer.sample_explorer import Rag_embedding, Transcriptome_embedding
from sample_explorer.utils import MsigDB_store, Sample_to_series_map
import pandas as pd
from tqdm import tqdm
import traceback
from typing import List
import pickle
import os

DATA_DIR = "data/"
TARGET_DIR = "tempfiles/transcriptome_seed_semantic_expand/"
GENESET_DIR = "gene_sets/"

os.makedirs(TARGET_DIR, exist_ok=True)

# Initialize embeddings and stores

rag_index = pd.read_pickle(DATA_DIR + "rag_index_v2.pkl")
with open(DATA_DIR + "rag_embedding_matrix_v2.pkl", "rb") as f:
    rag_embedding_matrix = pickle.load(f)
rag = Rag_embedding(rag_index, rag_embedding_matrix)


transcription_index = pd.read_pickle(DATA_DIR + "transcription_index_v2.p")
with open(DATA_DIR + "transcription_embedding_matrix_v2.pkl", "rb") as f:
    transcription_embedding_matrix = pickle.load(f)
exp = Transcriptome_embedding(transcription_index, transcription_embedding_matrix)

msigdb_store = MsigDB_store(DATA_DIR + "msigdb.v2023.2.Hs.symbols.gmt", DATA_DIR + "test_set_misgdb.csv")

y = Sample_to_series_map(DATA_DIR + "whole_metadata_human.p", rag_index, rag_embedding_matrix)


def load_test_genesets() -> List[str]:
    test_gs = pd.read_csv(f"{GENESET_DIR}test_gene_sets.txt", header=None)[0].to_list()
    return test_gs

# Get gene lists
gene_lists_collection = msigdb_store.list_genesets()
filtered_list = load_test_genesets()

# Load enrichment data
comb = pd.read_pickle(DATA_DIR + "all_enrichment.p")
ind_col = len(comb.columns)

# Initialize accumulator list
accum_list = []

def process_signature(geneset_name, comb, msigdb_store, y, rag):
    res_dict = {}
    testset = msigdb_store.get_gene_set_by_name(geneset_name)
    top1000 = comb[geneset_name].sort_values().head(1000).index.to_list()
    query = testset["long_title"]
    
    filesave = f"{TARGET_DIR}{geneset_name}.csv"
    
    series_of_interest = y.return_series_for_samples(top1000)
    additional_series = pd.concat([
        rag.get_closest_semantic_studies(series_id, k=6) 
        for series_id in series_of_interest["series_id"].unique()
    ])
    additional_series.drop_duplicates().to_csv(filesave)
    
    all_series = additional_series["gse_id_text"].drop_duplicates()
    #all_series_df = pd.DataFrame(all_series, columns=["gse_id"])
    #all_series_df.to_csv(f"{TARGET_DIR}{geneset_name}_all_series.csv", index=False)

    print(all_series)
    
    res_dict["signature"] = geneset_name
    res_dict["sample_explorer_sem"] = rag.get_averages_between_queries(query, all_series.to_list())
    
    return res_dict

for i in tqdm(filtered_list, desc="Processing signatures", unit="signature"):
    try:
        res_dict = process_signature(i, comb, msigdb_store, y, rag)
        accum_list.append(res_dict)
    except Exception as e:
        print(f"Failure processing batch {i}: {str(e)}")
        print(traceback.format_exc())

pd.DataFrame(accum_list).to_csv("tempfiles/cosine_similarity_tran_and_sem_queries.csv")

