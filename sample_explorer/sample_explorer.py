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
import sys
warnings.filterwarnings('ignore')
 
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

class Query_DB:
    logger = logging.getLogger(__name__)
    
    def __init__(self, semantic_vector_store, transcriptomic_vector_store, h5file=None):
        self.logger = logging.getLogger(__name__)
        H = logging.StreamHandler(sys.stdout)
        H.setLevel(logging.INFO)
        H.setFormatter(
            logging.Formatter(
                fmt="[%(asctime)s] %(levelname)s: %(message)s",
                datefmt="%d/%m/%Y ( %H:%M:%S )"
            ))
        self.logger.addHandler(H)

        self.logger.info("Loading transcriptomic vector store...")
        trans_obj = ad.read_h5ad(transcriptomic_vector_store, backed="r")
        self.logger.info("Transcriptomic vector store loaded.")
        
        self.logger.info("Loading semantic vector store...")
        sem_obj = ad.read_h5ad(semantic_vector_store, backed="r")
        self.logger.info("Semantic vector store loaded.")
        
        self.transcriptome_embedding = Transcriptome_embedding(trans_obj.obs, trans_obj.obsm["embedding"])
        self.logger.info("Transcriptome embedding initialized.")
        self.rag_embedding = Rag_embedding(sem_obj.obs, sem_obj.X)
        self.logger.info("RAG embedding initialized.")
        self.transcriptome_enrichment = Transcriptome_enrichment(trans_obj)        
        self.logger.info("Transcriptome object initialized.")
        
        if h5file is not None:
            self.logger.info("Loading ARCHS4 database object...")
            self.RNASeqAnalysis = RNASeqAnalysis(h5file)
            self.h5file = h5file    
            self.logger.info("ARCHS4 database object loaded.")
            
            self.logger.info("Creating series to sample mapping...")
            self.metafile = self.process_metadata(h5file)
            self.logger.info("DONE.")
        else:
            self.logger.info("ARCHS4 database object not loaded.")
            self.RNASeqAnalysis = None
        
        self.logger.info("Query_DB object initialized.")

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

    def retrieve_sample_metadata_from_h5(self, h5_path, sample_list):
        """
        Retrieves sample metadata from an HDF5 file.

        Args:
            h5_path (str): The path to the HDF5 file.
            sample_list (list): A list of sample names.

        Returns:
            dict: A dictionary containing the sample metadata.

        """
        sample_meta = a4.meta.samples(h5_path, sample_list)
        return sample_meta

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
        user_batch_size = 500
        output_df1_list = []
        self.logger.info("Starting transcriptome search...")
        for i in tqdm(range(0, len(my_list), user_batch_size)):
            try:
                samps = my_list[i:i + user_batch_size]
                df1 = self.transcriptome_enrichment.run_decouplr_on_memmaped_adata_with_samples(geneset, samps)
                # Store the output from each iteration
                output_df1_list.append(df1)
            except Exception as e:
                self.logger.error(f"Failure in batch {i}: {str(e)}")
                self.logger.error(traceback.format_exc())
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
        self.logger.info("Starting semantic search...")
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
        self.logger.info("Performing transcriptome expansion...")
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
        self.logger.info("Performing semantic expansion...")
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
        self.logger.info("Starting transcriptome search with semantic expansion...")
        self.logger.warning("This search strategy has long runtimes.")
        self.logger.info("search: " + str(search)),
        self.logger.info("expand: " + str(expand))
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
        self.logger.info("Starting transcriptome search with transcriptome expansion...")
        self.logger.warning("This search strategy has long runtimes.")
        self.logger.info("search: " + str(search))
        self.logger.info("expand: " + str(expand))  
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
        self.logger.info("Starting semantic search with semantic expansion...")
        self.logger.info("search: " + str(search))  
        self.logger.info("expand: " + str(expand))
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
        self.logger.info("Starting semantic search with transcriptome expansion...")
        self.logger.info("search: " + str(search))
        self.logger.info("expand: " + str(expand))
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
        if geneset is None and text_query is None:
            raise ValueError("At least one of geneset or text_query should not be None.")
        
        if (n_seed is None and n_expansion is not None) or (n_seed is not None and n_expansion is None):
            raise ValueError("Both n_seed and n_expansion must be either None or filled.")
        
        if n_seed == 0:
            raise ValueError("n_seed should not be zero.")

        if geneset is not None:
            if len(geneset) < 5:
                self.logger.warning("Gene set should have at least 5 genes.")
        
            if (self.transcriptome_enrichment.memmap_adata.var["gene"].isin(geneset)).sum() == 0:
                raise ValueError("No genes from the gene set found in the transcriptome.")
            
            if (self.transcriptome_enrichment.memmap_adata.var["gene"].isin(geneset)).sum() < len(geneset):
                self.logger.warning("Some genes from the gene set were not found in the transcriptome.")
                missing_genes = set(geneset) - set(self.transcriptome_enrichment.memmap_adata.var["gene"])
                missing_genes_percentage = (len(missing_genes) / len(geneset)) * 100
                self.logger.warning(f"Missing genes from geneset in transcriptome: {missing_genes} ({missing_genes_percentage:.2f}% of geneset)")

        results_object = Results(None, None, None) # new results object

        if text_query is None:
            search = "transcriptome"
            self.logger.info("Search type set to transcriptome by default.")
            if geneset is None or len(geneset) == 0:
                raise ValueError("Gene set should not be empty.")

        if geneset is None:
            search = "semantic"
            self.logger.info("Search type set to semantic by default.")
            if text_query is None or len(text_query) == 0:
                raise ValueError("Text query should not be empty.")
            if perform_enrichment:
                self.logger.error("Enrichment not performed. No gene set provided.")
                raise ValueError("Please specify a gene set.")

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

        self.logger.info("Completed search step.")


        if perform_enrichment and self.RNASeqAnalysis is not None:
            if results_object.expansion_studies is not None:
                self.logger.info("Performing enrichment on expansion studies.")
                samps = self.metafile[self.metafile["series_id"].isin(results_object.expansion_studies["gse_id"])]
                enrichment_df = self.RNASeqAnalysis.perform_enrichment_on_samples_batched(samps.index, geneset)
                res_df = pd.concat([enrichment_df, samps], axis = 1)
                results_object.samples = res_df
                self.logger.info("Completed enrichment on studies in expansion step.")
            else:
                samps = self.metafile[self.metafile["series_id"].isin(results_object.seed_studies["gse_id"])]
                self.logger.info("No expansion studies found. Performing enrichment on seed studies.")
                enrichment_df = self.RNASeqAnalysis.perform_enrichment_on_samples_batched(samps.index, geneset)
                res_df = pd.concat([enrichment_df, samps], axis = 1)
                results_object.samples = res_df
                self.logger.info("Completed enrichment on seed studies.")
        elif perform_enrichment is False and self.RNASeqAnalysis is not None:
            if results_object.expansion_studies is not None:
                samps = self.metafile[self.metafile["series_id"].isin(results_object.expansion_studies["gse_id"])]
                res_df = self.retrieve_sample_metadata_from_h5(self.h5file, samps["samples"])
                results_object.samples = res_df
                self.logger.info("Retrieving samples from expansion step. Enrichment not performed.")
            else:
                samps = self.metafile[self.metafile["series_id"].isin(results_object.seed_studies["gse_id"])]
                res_df = self.retrieve_sample_metadata_from_h5(self.h5file, samps["samples"])
                results_object.samples = res_df
                self.logger.info("Retrieving samples from search step. Enrichment not performed.")
        else:
            self.logger.info("Enrichment not performed.")
            results_object.samples = None
            
        #additional_series = additional_series.drop(["similarity_score"], axis=1).reset_index(drop=True).drop_duplicates() 
        #seed_series = seed_series.drop_duplicates()
        self.logger.info("query result retrieved.")
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
    



        

