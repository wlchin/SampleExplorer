from sample_explorer.sample_explorer import Transcriptome_embedding, Rag_embedding, Transcriptome_enrichment
import pandas as pd
from tqdm import tqdm
import traceback
from sample_explorer.utils import Sample_to_series_map
import pandas as pd
from sample_explorer.utils import MsigDB_store
import os
import pickle
from typing import List

DATA_DIR = "data/"
TARGET_DIR = "tempfiles/transcriptome_only/"
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

comb = pd.read_pickle("data/all_enrichment.p")

def load_test_genesets() -> List[str]:
    test_gs = pd.read_csv(f"{GENESET_DIR}test_gene_sets.txt", header=None)[0].to_list()
    return test_gs


filtered_list = load_test_genesets()
accum_list = []
for i in tqdm(filtered_list, desc="Processing signatures", unit="signature"):

    try:
        res_dict = {}
        testset = msigdb_store.get_gene_set_by_name(i)
        top1000 = comb[i].sort_values().head(1000).index.to_list()
        query = testset["long_title"]

        filesave = TARGET_DIR + i + ".csv"
        comb[i].sort_values().head(1000).to_csv(filesave)
        #res_dict["arch_samp"] = ind_arch
        res_dict["sample_explorer_sem"] = y.calculate_similarites_from_samples(query, top1000)
        res_dict["geneset"] = i
        accum_list.append(res_dict)

    except:
        print("failure", "batch: ", str(i))
        print(traceback.print_exc())


pd.DataFrame(accum_list).to_csv("tempfiles/cosine_similarity_tran_only_queries.csv")