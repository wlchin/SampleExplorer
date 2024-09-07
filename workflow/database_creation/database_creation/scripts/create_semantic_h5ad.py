import pickle
import anndata as ad
import pandas as pd

# Specify the path to your pickle file
pickle_file_path = "results/rag_embedding_matrix_v2.pkl"

# Open the pickle file in binary mode (rb) to read
with open(pickle_file_path, "rb") as f:
    # Load the data from the pickle file
    data = pickle.load(f)

dfs = pd.read_pickle("results/rag_index_v2.pkl").reset_index(drop = True)

adata = ad.AnnData(data)

adata.obs = dfs

adata.write("results/semantic_db.h5ad")