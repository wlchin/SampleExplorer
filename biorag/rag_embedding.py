import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from sentence_transformers import SentenceTransformer
import numpy as np


class Rag_embedding:
    def __init__(self, rag_index, rag_embedding_matrix):
        print("loading transformer")
        self.model = SentenceTransformer("paraphrase-MiniLM-L6-v2")
        self.rag_embedding_index = rag_index
        self.rag_embedding_matrix = rag_embedding_matrix

    def get_gse_from_rag_index(self, ind):
        """helper"""
        gse_id = self.rag_embedding_index.iloc[ind, :].to_list()[0]
        return gse_id

    def get_gse_text_from_rag_index(self, ind):
        gse_id_text = self.rag_embedding_index.iloc[ind, :].to_list()[1]
        return gse_id_text

    def get_text_linked_to_gse(self, gse_id_query):
        gse_id_text = self.rag_embedding_index[
            self.rag_embedding_index["series_id"].isin([gse_id_query])
        ]["rag_text"].to_list()[0]
        return gse_id_text

    def get_closest_semantic_studies(self, query_gse_id, k=5):
        """query against gse_id for expansion"""
        relevant_embedding_indices = np.where(
            self.rag_embedding_index["series_id"].isin([query_gse_id])
        )
        numpy_array = relevant_embedding_indices[0]
        flattened_array = numpy_array.flatten()
        query_embeddings = self.rag_embedding_matrix[flattened_array,]

        cosine_similarities = cosine_similarity(
            self.rag_embedding_matrix, query_embeddings
        )
        average_cosine_similarities = np.mean(cosine_similarities, axis=1)

        sorted_indices = np.argsort(average_cosine_similarities)[::-1]

        top_k = k
        results_data = []

        for i in range(top_k):
            index = sorted_indices[i]
            similarity = average_cosine_similarities[index]
            results_data.append({"Index": index, "similarity_score": similarity})

        # Create a Pandas DataFrame from the results_data
        results_df = pd.DataFrame(results_data)
        results_df["gse_id"] = [
            self.get_gse_from_rag_index(i) for i in results_df.Index
        ]
        results_df["gse_id_text"] = [
            self.get_gse_text_from_rag_index(i) for i in results_df.Index
        ]
        # results_df["similarity_source"] = "transcriptome"
        return results_df
    
    def query_rag(self, query, k=5):
        """query against text"""
        query_embedding = self.model.encode(query).reshape(1, -1)
        cosine_similarities = cosine_similarity(
            self.rag_embedding_matrix, query_embedding
        )
        cosine_distances = 1 - cosine_similarities
        sorted_indices = np.argsort(cosine_distances.flatten())

        top_k = k
        results_data = []

        for i in range(top_k):
            index = sorted_indices[i]
            distance = cosine_distances[index]
            results_data.append({"Index": index, "similarity_score": distance[0]})

        # Create a Pandas DataFrame from the results_data
        results_df = pd.DataFrame(results_data)
        results_df["gse_id"] = [self.get_gse_from_rag_index(i) for i in results_df.Index]
        results_df["gse_id_text"] = [
            self.get_gse_text_from_rag_index(i) for i in results_df.Index
        ]
        results_df["similarity_source"] = "semantic"

        return results_df
    
    def get_averages_between_queries(self, query, list_of_queries):
        query_embedding = self.model.encode([query])[0]
        list_embeddings = self.model.encode(list_of_queries)
        
        query_embedding = query_embedding.reshape(1, -1)
        list_embeddings = list_embeddings.reshape(len(list_of_queries), -1)

        # Calculate cosine similarity using scikit-learn
        cosine_scores = cosine_similarity(query_embedding, list_embeddings)

        # Calculate average similarity
        average_similarity = np.mean(cosine_scores)
        return average_similarity