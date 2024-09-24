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
TARGET_DIR = "tempfiles/transcriptome_seed_transcriptome_expand/"
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

# files_os = os.listdir("enrichment")
# selected_strings = [string for string in files_os if 'pval' in string]
# selpath = ["enrichment/" + i for i in selected_strings]
# dfs = [pd.read_pickle(i) for i in selpath]
comb = pd.read_pickle("data/all_enrichment.p")

def get_series_of_relevance_from_series(series_of_interest, exp, k = 5):
    series_list = []
    target_series = series_of_interest["series_id"].drop_duplicates()
    for i in tqdm(target_series):
        res_temp = exp.get_closest_transcriptional_studies(i, k = k)
        series_list.append(res_temp)
    additional_series = pd.concat(series_list).drop_duplicates()
    return additional_series

def save_file(res, directory, filename):
    file_path = os.path.join(directory, filename)
    res.to_csv(file_path)

accum_list = []

filtered_list = load_test_genesets()

for i in tqdm(filtered_list, desc="Processing signatures", unit="signature"):
    try:
        res_dict = {}
        testset = msigdb_store.get_gene_set_by_name(i)
        top1000 = comb[i].sort_values().head(1000).index.to_list()
        series_of_interest = y.return_series_for_samples(top1000)
        query = testset["long_title"]
        additional_series = get_series_of_relevance_from_series(series_of_interest, exp) # change
        dest_folder = TARGET_DIR
        filename = i + ".csv"
        save_file(additional_series, dest_folder, filename)
        all_series = additional_series["gse_id"].drop_duplicates().tolist()
        res_dict["signature"] = i
        res_dict["sample_explorer_sem"] = y.calculate_similarites_from_series(query, all_series)
        accum_list.append(res_dict)
    except:
        print("failure", "batch: ", str(i))
        print(traceback.print_exc())

pd.DataFrame(accum_list).to_csv("tempfiles/cosine_similarity_tran_and_trans_queries.csv")

