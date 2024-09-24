import numpy as np
import pandas as pd
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
from sample_explorer.utils import MsigDB_store, Sample_to_series_map
import os


metadata = Sample_to_series_map("data/whole_metadata_human.p", 
                                "data/important_data_for_LLM/rag_index_v2.pkl", 
                                "data/important_data_for_LLM/rag_embedding_matrix_v2.pkl")
rag_index = pd.read_pickle("data/important_data_for_LLM/rag_index_v2.pkl")
msigdb_store = MsigDB_store("data/msigdb.v2023.2.Hs.symbols.gmt", "data/test_set_misgdb.csv")

def calculate_cosine_similarity(input_string, list_of_strings):
    # Load a pre-trained model from sentence-transformers
    model = SentenceTransformer('paraphrase-MiniLM-L6-v2')

    # Calculate the embeddings for the input string and the list of strings
    input_embedding = model.encode(input_string, convert_to_tensor=True).reshape(1, -1)
    list_embeddings = model.encode(list_of_strings, convert_to_tensor=True)

    # Convert the list of embeddings to a 2D NumPy array
    list_embeddings = list_embeddings.cpu().detach().numpy()

    # Calculate cosine similarity
    similarity_scores = cosine_similarity(input_embedding, list_embeddings)

    return np.mean(similarity_scores.flatten())




store_df = []
BASE_DIR = "tempfiles/arch_enrichment/"
file_list = os.listdir(BASE_DIR)

for i in file_list: 
    arch_gs_query = i.replace(".csv", "")
    gene_set_dict = msigdb_store.get_gene_set_by_name(arch_gs_query)
    long_query = gene_set_dict["long_title"]
    archsamp = pd.read_csv(BASE_DIR + i, index_col=0)
    samps = archsamp.index.to_list()
    df_series_samples = metadata.return_series_for_samples(samps)

    to_consider_ARCH = rag_index[rag_index["series_id"].isin(df_series_samples["series_id"])]
    arch_metric = calculate_cosine_similarity(long_query, to_consider_ARCH["rag_text"].to_list())
    geneset_dict = {"geneset": arch_gs_query, "arch_metric": arch_metric}
    store_df.append(geneset_dict)


df = pd.DataFrame(store_df)
df.to_csv("tempfiles/cosine_similarity_arch_queries.csv")

