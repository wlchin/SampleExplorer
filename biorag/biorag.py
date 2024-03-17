import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from sentence_transformers import SentenceTransformer
import anndata as ad
import decoupler as dc
import numpy as np
from scipy.spatial.distance import pdist
import statsmodels.stats.multitest as smm
import pickle
import archs4py as a4
import os
import warnings

os.environ["TOKENIZERS_PARALLELISM"] = "false"

warnings.filterwarnings('ignore')
 
class Rag_embedding_adaptor:
    def __init__(self, h5ad_path):
        self.rag_index = pd.read_pickle(rag_index_path)
        self.rag_embedding = pd.read_pickle(rag_embedding_path)
        
    def get_text_linked_to_gse(self, gse_id):
        return self.rag_index.loc[gse_id]["text"]
    
    def get_averages_between_queries(self, query, experimental_text):
        query_embedding = self.rag_embedding.encode(query)
        experimental_text_embeddings = self.rag_embedding.encode(experimental_text)
        return np.mean(cosine_similarity(query_embedding, experimental_text_embeddings))
    
class Transcriptome_embedding_adaptor:
    def __init__(self, h5ad_path):
        pass

class BioRAG:
    def __init__(self, h5_path, rag_index_path, rag_embedding_path):
        self.geo_metadata = self.process_metadata(h5_path)
        self.rag_embedding_object = Rag_embedding_adaptor(rag_index_path, rag_embedding_path) 

        
