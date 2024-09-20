import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from sentence_transformers import SentenceTransformer
import numpy as np
import logging 
import sys

class Rag_embedding:
    def __init__(self, rag_index, rag_embedding_matrix):
        self.logger = logging.getLogger(__name__)
        H = logging.StreamHandler(sys.stdout)
        H.setLevel(logging.INFO)
        H.setFormatter(
            logging.Formatter(
                fmt="[%(asctime)s] %(levelname)s: %(message)s",
                datefmt="%d/%m/%Y ( %H:%M:%S )"
            ))
        self.logger.addHandler(H)


        self.model = SentenceTransformer("paraphrase-MiniLM-L6-v2")
        self.logger.info("Loaded SentenceTransformer model.")
        self.rag_embedding_index = rag_index
        self.rag_embedding_matrix = rag_embedding_matrix

    def get_gse_from_rag_index(self, ind):
        """Get the GSE ID from the RAG embedding index.

        Args:
            ind (int): The index of the GSE ID in the RAG embedding index.

        Returns:
            str: The GSE ID corresponding to the given index.
        """
        gse_id = self.rag_embedding_index.iloc[ind, :].to_list()[0]
        return gse_id

    def get_gse_text_from_rag_index(self, ind):
        """
        Retrieves the GSE text from the RAG index.

        Parameters:
        - ind (int): The index of the GSE text to retrieve.

        Returns:
        - gse_id_text (str): The GSE text corresponding to the given index.
        """
        gse_id_text = self.rag_embedding_index.iloc[ind, :].to_list()[1]
        return gse_id_text

    def get_text_linked_to_gse(self, gse_id_query):
        """
        Retrieves the text linked to a given GSE ID.

        Parameters:
        - gse_id_query (str): The GSE ID to query.

        Returns:
        - str: The text linked to the specified GSE ID.
        """
        gse_id_text = self.rag_embedding_index[
            self.rag_embedding_index["series_id"].isin([gse_id_query])
        ]["rag_text"].to_list()[0]
        return gse_id_text

    def get_closest_semantic_studies(self, query_gse_id, k=5):
        """
        Retrieves the closest semantic studies based on a given query GSE ID.

        Parameters:
            query_gse_id (str): The GSE ID to query against for expansion.
            k (int): The number of closest studies to retrieve. Default is 5.

        Returns:
            pandas.DataFrame: A DataFrame containing the closest studies along with their similarity scores.

        """
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
        """
        Query the RAG (Retrieval-Augmented Generation) model with a given query.

        Args:
            query (str): The query text.
            k (int, optional): The number of top results to retrieve. Defaults to 5.

        Returns:
            pandas.DataFrame: A DataFrame containing the top k results with columns 'similarity_score', 'gse_id', 'gse_id_text', and 'similarity_source'.
        """
        query_embedding = self.model.encode(query).reshape(1, -1)
        cosine_similarities = cosine_similarity(
            self.rag_embedding_matrix, query_embedding
        )
        cosine_distances = 1 - cosine_similarities
        sorted_indices = np.argsort(cosine_distances.flatten())

        if len(sorted_indices) < k:
            k = len(sorted_indices)

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

        return results_df.drop("Index", axis=1)
    
    def get_averages_between_queries(self, query, list_of_queries):
        """
        Calculates the average cosine similarity between a query and a list of queries.

        Parameters:
        query (str): The query string.
        list_of_queries (list): A list of query strings.

        Returns:
        float: The average cosine similarity between the query and the list of queries.
        """
        query_embedding = self.model.encode([query])[0]
        list_embeddings = self.model.encode(list_of_queries)
        
        query_embedding = query_embedding.reshape(1, -1)
        list_embeddings = list_embeddings.reshape(len(list_of_queries), -1)

        # Calculate cosine similarity using scikit-learn
        cosine_scores = cosine_similarity(query_embedding, list_embeddings)

        # Calculate average similarity
        average_similarity = np.mean(cosine_scores)
        return average_similarity