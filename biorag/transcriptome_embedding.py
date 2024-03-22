import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
from scipy.spatial.distance import pdist


class Transcriptome_embedding:
    def __init__(self, embedding_index, embedding_matrix):
        self.embeddings_index = embedding_index
        self.embedding_matrix = embedding_matrix

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