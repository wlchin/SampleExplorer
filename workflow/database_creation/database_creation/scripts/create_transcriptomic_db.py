import anndata as ad
import numpy as np
import pandas as pd
import pickle


rows = 287553
cols = 67186

memmap_filename = 'results/study_counts.dat'
reloaded_memmap_matrix = np.memmap(memmap_filename, dtype='float16', mode='r', shape=(rows, cols))

test = np.copy(reloaded_memmap_matrix)
adata = ad.AnnData(test)
obs_df = pd.read_pickle("results/transcription_index_v2.pkl")
var_df = pd.read_pickle("results/genes_v2.p")
var_df.columns = ["gene"]
adata.obs = obs_df
adata.var = var_df

with open("results/transcription_embedding_matrix_v2.pkl", 'rb') as file:
    embedding_matrix = pickle.load(file) 

adata.obsm["embedding"] = embedding_matrix
adata.write('results/transcriptomic_db.h5ad')

