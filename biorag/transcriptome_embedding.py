import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from sentence_transformers import SentenceTransformer
import anndata as ad
import decoupler as dc
import numpy as np
from scipy.spatial.distance import pdist
import statsmodels.stats.multitest as smm
import pickle

class Transcriptome_embedding:
    def __init__(self, embedding_index_path, embedding_matrix_path):
        self.embeddings_index = pd.read_pickle(embedding_index_path)
        with open(embedding_matrix_path, 'rb') as file:
            self.embedding_matrix = pickle.load(file) 

    def search_relevant_embeddings_by_sample(self, samples):
        searching = self.embedding_matrix[self.embeddings_index.index.isin(samples),:]
        return searching

    def search_relevant_embeddings_by_series(self, series):
        samples = self.embeddings_index[self.embeddings_index["series_id"].isin(series)].index
        searching = self.embedding_matrix[self.embeddings_index.index.isin(samples),:]
        return searching
    
    def obtain_samples_from_series(self, series):
        samples = self.embeddings_index[self.embeddings_index["series_id"].isin(series)].index
        return samples
    
    def mean_interpoint_distance(self, data, p = 1):
    # Calculate pairwise distances between points using Minkowski distance
        distances = pdist(data, metric='minkowski', p=p)

        # Calculate the mean distance
        mean_distance = np.mean(distances)
        scientific_notation = "{:.2e}".format(mean_distance)

        return scientific_notation
    
    def get_gse_from_transcriptome_index(self, ind):
        gse_id = self.embeddings_index.iloc[ind,:].to_list()[0]
        return gse_id
    
    def calculate_interpoint_distance_for_samples(self, samples, p):
        searching = self.search_relevant_embeddings_by_sample(samples)
        dist = self.mean_interpoint_distance(searching, p)
        return dist

    def calculate_interpoint_distance_for_series(self, series, p):
        searching = self.search_relevant_embeddings_by_series(series)
        dist = self.mean_interpoint_distance(searching, p)
        return dist
    
    def get_closest_transcriptional_studies(self, query_gse_id, k = 5):
        relevant_embedding_indices = np.where(self.embeddings_index["series_id"].isin([query_gse_id]))
        numpy_array = relevant_embedding_indices[0]
        flattened_array = numpy_array.flatten()
        query_embeddings = self.embedding_matrix[flattened_array,]
        
        cosine_similarities = cosine_similarity(self.embedding_matrix, query_embeddings)
        average_cosine_similarities = np.mean(cosine_similarities, axis=1)
        
        sorted_indices = np.argsort(average_cosine_similarities)[::-1]
        
        top_k = k
        results_data = []
    
        for i in range(top_k):
            index = sorted_indices[i]
            similarity = average_cosine_similarities[index]
            results_data.append({'Index': index, 'similarity_score': similarity})
    
        # Create a Pandas DataFrame from the results_data
        results_df = pd.DataFrame(results_data)
        results_df["gse_id"] = [self.get_gse_from_transcriptome_index(i) for i in results_df.Index]
        return results_df

    def get_closest_transcriptional_studies_by_sample(self, query_sample, k = 5):
        embedding = self.search_relevant_embeddings_by_sample(query_sample)
        query_vector = embedding.reshape(1, -1)
        cosine_similarities = cosine_similarity(query_vector, self.embedding_matrix)
        closest_indices = np.argsort(cosine_similarities[0])[-k:][::-1]
        closest_samples = self.embeddings_index.iloc[closest_indices,:]
        return closest_samples
# Calculate cosine similarity

class Transcriptome_enrichment:
    def __init__(self, counts_path, gene_path, transcription_index_path):
        self.memmap_adata = self.create_memmap_adata(counts_path, gene_path, transcription_index_path)
        
    def create_memmap_adata(self, counts_path, gene_path, transcription_index_path):
        rows = 167649
        cols = 39376
        memmap_filename = counts_path
        reloaded_memmap_matrix = np.memmap(memmap_filename, dtype='float32', mode='r', shape=(rows, cols))
        adata = ad.AnnData(reloaded_memmap_matrix)
        obs_df = pd.read_pickle(transcription_index_path)
        var_df = pd.read_pickle(gene_path)
        adata.obs = obs_df
        adata.var = var_df
        return adata
    
    def list_to_dc_geneset(self, list_test_geneset):
        df = pd.DataFrame(list_test_geneset)
        df.columns = ["genesymbol"]
        df["geneset"] = "query_gene_set"
        return df
    
    def list_to_dc_geneset_dictionary(self, list_test_geneset_dict):
        geneset_df_list = []
        for i,j in list_test_geneset_dict.items():
            df = pd.DataFrame(j)
            df.columns = ["genesymbol"]
            df["geneset"] = i
            geneset_df_list.append(df)
        return pd.concat(geneset_df_list)
    
  
    def run_decouplr_on_memmaped_adata(self, genelist, series_list, nstudies = 10):
        
        if isinstance(genelist, dict):
            geneset = self.list_to_dc_geneset_dictionary(genelist)
            #print(geneset)

        elif isinstance(genelist, list):
            geneset = self.list_to_dc_geneset(genelist)
            #print(geneset)
        
        else:
            "incorrect format"

        adata = self.memmap_adata
        query_adata = adata[adata.obs["series_id"].isin(series_list)]
        
        dc.run_ora(
            mat=query_adata,
            net=geneset,
            source='geneset',
            target='genesymbol',
            verbose=True, use_raw=False
        )
        if query_adata.obsm["ora_pvals"].shape[1] > 1:
            return query_adata.obsm["ora_pvals"], query_adata.obsm["ora_estimate"]
        else:
            res = pd.concat([query_adata.obsm["ora_pvals"], query_adata.obsm["ora_estimate"]], axis = 1)
            res.columns = ["pvals", "estimate"]
            _, p_corrected, _, _ = smm.multipletests(res['pvals'], method='fdr_bh')
            res['p_corrected'] = p_corrected
            #studies = res.sort_values("estimate", ascending = False).head(nstudies).index.str.replace("series_id.", "")
            return res
        
    def run_decouplr_on_memmaped_adata_with_samples(self, genelist, samples_list, nstudies = 10):

        if isinstance(genelist, dict):
            geneset = self.list_to_dc_geneset_dictionary(genelist)
            #print(geneset)

        elif isinstance(genelist, list):
            geneset = self.list_to_dc_geneset(genelist)
            #print(geneset)

        else:
            "incorrect format"

        adata = self.memmap_adata
        query_adata = adata[adata.obs.index.isin(samples_list)]

        dc.run_ora(
            mat=query_adata,
            net=geneset,
            source='geneset',
            target='genesymbol',
            verbose=True, use_raw=False
        )
        if query_adata.obsm["ora_pvals"].shape[1] > 1:
            return query_adata.obsm["ora_pvals"], query_adata.obsm["ora_estimate"]
        else:
            res = pd.concat([query_adata.obsm["ora_pvals"], query_adata.obsm["ora_estimate"]], axis = 1)
            res.columns = ["pvals", "estimate"]
            _, p_corrected, _, _ = smm.multipletests(res['pvals'], method='fdr_bh')
            res['p_corrected'] = p_corrected
            return res
    
        
