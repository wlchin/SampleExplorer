import os
import pickle
import traceback
from typing import List
import pandas as pd

from sample_explorer.sample_explorer import Rag_embedding
from sample_explorer.utils import MsigDB_store

# Constants
DATA_DIR = "data"
TEMP_DIR = "tempfiles/sem_query_rag_df"
GENESET_DIR = "gene_sets"

def load_test_genesets() -> List[str]:
    test_gs = pd.read_csv(f"{GENESET_DIR}/test_gene_sets.txt", header=None)[0].to_list()
    return test_gs

def check_folder_existence():
    os.makedirs(TEMP_DIR, exist_ok=True)

def load_rag_data():
    rag_index = pd.read_pickle(f"{DATA_DIR}/important_data_for_LLM/rag_index_v2.pkl")
    with open(f"{DATA_DIR}/important_data_for_LLM/rag_embedding_matrix_v2.pkl", "rb") as f:
        rag_embedding_matrix = pickle.load(f)
    return Rag_embedding(rag_index, rag_embedding_matrix)

def process_gene_set(rag: Rag_embedding, msigdb_store: MsigDB_store, gene_set_name: str):
    check_folder_existence()
    try:
        testset = msigdb_store.get_gene_set_by_name(gene_set_name)
        test_text = testset["long_title"]
        df_res = rag.query_rag(test_text, k=50)
        savefile = os.path.join(TEMP_DIR, f"{gene_set_name}.csv")
        df_res.to_csv(savefile)
    except Exception as e:
        print(f"Failure: {gene_set_name}")
        traceback.print_exc()

def main():
    # Initialize data
    test_gene_sets = load_test_genesets()
    rag = load_rag_data()
    msigdb_store = MsigDB_store(f"{DATA_DIR}/msigdb.v2023.2.Hs.symbols.gmt", f"{DATA_DIR}/test_set_misgdb.csv")

    for i in test_gene_sets:
        process_gene_set(rag, msigdb_store, i)


if __name__ == "__main__":
    main()
