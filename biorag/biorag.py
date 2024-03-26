import pandas as pd
from tqdm import tqdm
import anndata as ad
import archs4py as a4
import os
import traceback
from .rnaseq_analysis import RNASeqAnalysis
from .rag_embedding import Rag_embedding
from .transcriptome_embedding import Transcriptome_embedding
from .enrichment import Transcriptome_enrichment

os.environ["TOKENIZERS_PARALLELISM"] = "false"
import warnings
import logging
import time
warnings.filterwarnings('ignore')
 

 
class Query_DB:
    def __init__(self, semantic_vector_store, transcriptomic_vector_store, h5file=None, log_level=logging.INFO):
        logging.basicConfig(filename='log.txt', level=log_level, format='%(asctime)s %(levelname)s %(message)s')
        logging.info("Loading transcriptomic vector store...")
        trans_obj = ad.read_h5ad(transcriptomic_vector_store, backed="r")
        logging.info("Transcriptomic vector store loaded.")
        
        logging.info("Loading semantic vector store...")
        sem_obj = ad.read_h5ad(semantic_vector_store, backed="r")
        logging.info("Semantic vector store loaded.")
        
        self.transcriptome_embedding = Transcriptome_embedding(trans_obj.obs, trans_obj.obsm["embedding"])
        logging.info("Transcriptome embedding initialized.")
        self.rag_embedding = Rag_embedding(sem_obj.obs, sem_obj.X)
        logging.info("RAG embedding initialized.")
        self.transcriptome_enrichment = Transcriptome_enrichment(trans_obj)        
        logging.info("Transcriptome object initialized.")
        
        if h5file is not None:
            logging.info("Loading RNASeqAnalysis object...")
            self.RNASeqAnalysis = RNASeqAnalysis(h5file)
            logging.info("RNASeqAnalysis object loaded.")
            
            logging.info("Processing metadata...")
            self.metafile = self.process_metadata(h5file)
            logging.info("Metadata processed.")
        else:
            logging.info("RNASeqAnalysis object not loaded.")
            self.RNASeqAnalysis = None
            #self.sample_to_series_map = Sample_to_series_map(h5file, self.rag_embedding)

    def process_metadata(self, h5_path):
        """
        Process metadata from an HDF5 file.

        Args:
            h5_path (str): The path to the HDF5 file.

        Returns:
            pandas.DataFrame: A DataFrame containing processed metadata with columns 'series_id' and 'samples'.
        """
        x = a4.meta.meta(h5_path, ".*", meta_fields=["series_id"])
        x["samples"] = x.index
        x = x[['series_id', 'samples']].drop_duplicates()
        return x

    def get_top_samples(self, df, n=1000):
        """
        Returns a subset of samples from the given DataFrame based on their ranking.

        Parameters:
        - df (pandas.DataFrame): The DataFrame containing the samples.
        - n (int): The number of top samples to retrieve. Default is 1000.

        Returns:
        - pandas.DataFrame: A subset of the original DataFrame containing the top samples.
        """
        df_meta = self.transcriptome_embedding.embeddings_index.copy()
        samps = df.iloc[:, 1].sort_values(ascending=False).head(n).index.to_list()
        series_of_interest = df_meta[df_meta.index.isin(samps)]
        return series_of_interest

    def transcriptome_search(self, geneset, nsamples=1000):
        """
        Perform a transcriptome search using a given geneset.

        Args:
            geneset (list): A list of genes to search for in the transcriptome.
            nsamples (int, optional): The number of top samples to retrieve. Defaults to 1000.

        Returns:
            tuple: A tuple containing two pandas DataFrames. The first DataFrame contains the top samples of interest,
                   and the second DataFrame contains the relevant series (studies) associated with the top samples.
        """
        my_list = self.transcriptome_enrichment.memmap_adata.obs.index
        user_batch_size = 250
        output_df1_list = []
        for i in tqdm(range(0, len(my_list), user_batch_size), desc="Processing Batches", unit="batch"):
            try:
                samps = my_list[i:i + user_batch_size]
                df1 = self.transcriptome_enrichment.run_decouplr_on_memmaped_adata_with_samples(geneset, samps)
                # Store the output from each iteration
                output_df1_list.append(df1)
            except Exception as e:
                logging.error(f"Failure in batch {i}: {str(e)}")
                logging.error(traceback.format_exc())
        final_df1 = pd.concat(output_df1_list)
        series_of_interest = self.get_top_samples(final_df1, n=nsamples)
        relevant_series = series_of_interest["series_id"]  # extract only studies (series)
        series_of_interest = series_of_interest.rename(columns={"series_id": "gse_id"})
        return series_of_interest, relevant_series
        
    
    def semantic_search(self, test_text, k=50):
        """
        Perform semantic search using the RAG embedding model.

        Parameters:
        - test_text (str): The text to search for.
        - k (int): The number of results to retrieve (default is 50).

        Returns:
        - df_res (DataFrame): The search results as a DataFrame.
        - df_res["gse_id"] (Series): The series containing the "gse_id" column from the search results.
        """
        df_res = self.rag_embedding.query_rag(test_text, k=k)
        return df_res, df_res["gse_id"]
 
    def get_transcriptome_series_of_relevance_from_series(self, series_of_interest, k=5):
        """
        Retrieves a series of transcriptome data that are relevant to the given series of interest.

        Parameters:
        - series_of_interest (list): A list of series of interest.
        - k (int): The number of closest transcriptional studies to retrieve for each series.

        Returns:
        - additional_series (DataFrame): A DataFrame containing the additional series of transcriptome data.
        """

        series_list = []
        target_series = series_of_interest
        for i in tqdm(target_series):
            res_temp = self.transcriptome_embedding.get_closest_transcriptional_studies(i, k=k)
            series_list.append(res_temp)
        additional_series = pd.concat(series_list)
        additional_series["gse_id_text"] = [self.rag_embedding.get_text_linked_to_gse(i) for i in additional_series["gse_id"]]
        return additional_series
    
    def get_semantic_series_of_relevance_from_series(self, series_of_interest, k=5):
        """
        Retrieves additional series of relevance based on the input series of interest.

        Parameters:
        - series_of_interest (list): A list of series of interest.
        - k (int): The number of closest semantic studies to retrieve for each series.

        Returns:
        - additional_series (pandas.DataFrame): A DataFrame containing the additional series of relevance.

        Usage:
        Call this method to retrieve additional series of relevance based on the input series of interest.
        The method uses the rag_embedding object to get the closest semantic studies for each series.
        The number of closest semantic studies to retrieve can be specified using the k parameter.
        """
        additional_series = pd.concat([self.rag_embedding.get_closest_semantic_studies(i, k=k) for i in series_of_interest])
        return additional_series
    

    def transcriptome_search_with_semantic_expansion(self, geneset_query, search=1000, expand=5):
        """
        Perform a transcriptome search with semantic expansion.

        Args:
            geneset_query (str): The query for the geneset.
            search (int, optional): The number of samples to search. Defaults to 1000.
            expand (int, optional): The number of additional series to expand the search. Defaults to 5.

        Returns:
            tuple or None: A tuple containing the additional series and the series dataframe if expand is not 0,
            otherwise None and the series dataframe.
        """
        if expand == 0:
            series_df, series_of_interest = self.transcriptome_search(geneset_query, search)
            return None, series_df
        else:
            series_df, series_of_interest = self.transcriptome_search(geneset=geneset_query, nsamples=search)
            additional_series = self.get_semantic_series_of_relevance_from_series(series_of_interest, expand)
            return additional_series, series_df

    def transcriptome_search_with_transcriptome_expansion(self, geneset_query, search=1000, expand=10):
        """
        Perform a transcriptome search with transcriptome expansion.

        Args:
            geneset_query (str): The geneset query to search for.
            search (int): The number of search results to retrieve. Default is 1000.
            expand (int): The number of additional series to expand the search. Default is 10.

        Returns:
            tuple or None: A tuple containing the additional series and the series dataframe if expand is not 0,
            otherwise None and the series dataframe.

        """
        if expand == 0:
            series_df, series_of_interest = self.transcriptome_search(geneset_query, search)
            return None, series_df
        else:
            series_df, series_of_interest = self.transcriptome_search(geneset_query, search)
            additional_series = self.get_transcriptome_series_of_relevance_from_series(series_of_interest, expand)
            return additional_series, series_df

    def semantic_search_with_semantic_expansion(self, text_query, search=50, expand=5):
        """
        Perform semantic search with semantic expansion.

        Args:
            text_query (str): The text query to search for.
            search (int): The number of search results to retrieve. Default is 50.
            expand (int): The number of additional series to expand the search results. Default is 5.

        Returns:
            tuple or None: A tuple containing the additional series and the search results dataframe.
                           If `expand` is 0, returns None as the additional series.
        """
        print(text_query)
        if expand == 0:
            series_df, series_of_interest = self.semantic_search(text_query, search)
            return None, series_df
        else:
            series_df, series_of_interest = self.semantic_search(text_query, search)
            additional_series = self.get_semantic_series_of_relevance_from_series(series_of_interest, expand)
            return additional_series, series_df

    def semantic_search_with_transcriptome_expansion(self, text_query, search=50, expand=10):
        """
        Performs semantic search with transcriptome expansion.

        Args:
            text_query (str): The text query to search for.
            search (int): The number of search results to retrieve.
            expand (int): The number of additional transcriptome series to expand the search with.

        Returns:
            tuple or None: A tuple containing the additional transcriptome series and the search results dataframe.
                           If `expand` is 0, returns `None` as the first element of the tuple.
        """
        if expand == 0:
            series_df, series_of_interest = self.semantic_search(text_query, search)
            return None, series_df
        else:
            series_df, series_of_interest = self.semantic_search(text_query, search)
            additional_series = self.get_transcriptome_series_of_relevance_from_series(series_of_interest, expand)
            return additional_series, series_df

    def search(self, geneset, text_query, search="semantic", expand="transcriptome", perform_enrichment=False, n_seed=None, n_expansion=None):
        """
        Perform a search using the specified parameters.

        Args:
            geneset (list): A list of genes to search for.
            text_query (str): The text query to search for.
            search (str, optional): The type of search to perform. Defaults to "semantic".
            expand (str, optional): The type of expansion to perform. Defaults to "transcriptome".
            perform_enrichment (bool, optional): Whether to perform enrichment analysis. Defaults to False.
            n_seed (int, optional): The number of seed studies to include. Defaults to None.
            n_expansion (int, optional): The number of expansion studies to include. Defaults to None.

        Raises:
            ValueError: If the gene set is empty or has less than 5 genes.
            ValueError: If either n_seed or n_expansion is provided without the other.

        Returns:
            Results: An object containing the search results.
        """
        if geneset is None or len(geneset) == 0:
            raise ValueError("Gene set should not be empty.")
        
        if geneset is not None and len(geneset) < 5:
            logging.warning("Gene set should have at least 5 genes.")
        
        if (n_seed is None and n_expansion is not None) or (n_seed is not None and n_expansion is None):
            raise ValueError("Both n_seed and n_expansion must be either None or filled.")
        
        # Rest of the code...
    def search(self, geneset, text_query, search="semantic", expand="transcriptome", perform_enrichment=False, n_seed=None, n_expansion=None):
        """this is the main function
        """
        if geneset is None or len(geneset) == 0:
            raise ValueError("Gene set should not be empty.")
        
        if geneset is not None and len(geneset) < 5:
            logging.warning("Gene set should have at least 5 genes.")
        
        if (n_seed is None and n_expansion is not None) or (n_seed is not None and n_expansion is None):
            raise ValueError("Both n_seed and n_expansion must be either None or filled.")
        
        # Rest of the code...   # Rest of the code...

        results_object = Results(None, None, None) # new results object

        if text_query is None:
            search = "transcriptome"

        if geneset is None:
            search = "semantic"

        if search == "semantic":
            if expand == "transcriptome":
                if n_seed is None and n_expansion is None:
                    additional_series, seed_series = self.semantic_search_with_transcriptome_expansion(text_query)
                else:
                    additional_series, seed_series = self.semantic_search_with_transcriptome_expansion(text_query, search = n_seed, expand = n_expansion)
                results_object.seed_studies = seed_series
                results_object.expansion_studies = additional_series
            elif expand == "semantic":
                if n_seed is None and n_expansion is None:
                    additional_series, seed_series = self.semantic_search_with_semantic_expansion(text_query)
                else:
                    additional_series, seed_series = self.semantic_search_with_semantic_expansion(text_query, search = n_seed, expand = n_expansion)
                results_object.seed_studies = seed_series
                results_object.expansion_studies = additional_series
        if search == "transcriptome":
            if expand == "transcriptome":
                if n_seed is None and n_expansion is None:
                    additional_series, seed_series = self.transcriptome_search_with_transcriptome_expansion(geneset)
                else:
                    additional_series, seed_series = self.transcriptome_search_with_transcriptome_expansion(geneset, search = n_seed, expand = n_expansion)
                results_object.seed_studies = seed_series
                results_object.expansion_studies = additional_series
            elif expand == "semantic":
                if n_seed is None and n_expansion is None:
                    additional_series, seed_series = self.transcriptome_search_with_semantic_expansion(geneset)
                else :
                    additional_series, seed_series = self.transcriptome_search_with_semantic_expansion(geneset, search = n_seed, expand = n_expansion)
                results_object.seed_studies = seed_series
                results_object.expansion_studies = additional_series


        if perform_enrichment and self.RNASeqAnalysis is not None:
            if results_object.expansion_studies is not None:
                samps = self.metafile[self.metafile["series_id"].isin(results_object.expansion_studies["gse_id"])]
                enrichment_df = self.RNASeqAnalysis.perform_enrichment_on_samples_batched(samps.index, geneset)
                res_df = pd.concat([enrichment_df, samps], axis = 1)
                results_object.samples = res_df
            else:
                samps = self.metafile[self.metafile["series_id"].isin(additional_series["gse_id"])]
                res_df = samps
                results_object.samples = res_df
        else:
            results_object.samples = None
            
        #additional_series = additional_series.drop(["similarity_score"], axis=1).reset_index(drop=True).drop_duplicates() 
        #seed_series = seed_series.drop_duplicates()

        return results_object

class Results:
    """
    Represents the results of a bioRAG analysis.

    Attributes:
        seed_studies (DataFrame): The seed studies used in the analysis.
        expansion_studies (DataFrame): The expansion studies used in the analysis.
        samples (DataFrame): The samples used in the analysis.
    """

    def __init__(self, seed_df, expansion_df, samples_df):
        self.seed_studies = seed_df
        self.expansion_studies = expansion_df
        self.samples = samples_df
    



        

