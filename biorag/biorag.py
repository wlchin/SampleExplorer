import pandas as pd
from tqdm import tqdm
import anndata as ad
import numpy as np
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
warnings.filterwarnings('ignore')
 
 
class Query_DB:
    def __init__(self, semantic_vector_store, transcriptomic_vector_store, h5file = None):
        trans_obj = ad.read_h5ad(transcriptomic_vector_store, backed = "r")
        sem_obj = ad.read_h5ad(semantic_vector_store, backed = "r")
        self.transcriptome_embedding = Transcriptome_embedding(trans_obj.obs, trans_obj.obsm["embedding"])
        self.rag_embedding = Rag_embedding(sem_obj.obs, sem_obj.X)
        self.transcriptome_enrichment = Transcriptome_enrichment(trans_obj)
        if h5file is not None:
            self.RNASeqAnalysis = RNASeqAnalysis(h5file)
            self.metafile = self.process_metadata(h5file)
        else:
            self.RNASeqAnalysis = None
            #self.sample_to_series_map = Sample_to_series_map(h5file, self.rag_embedding)
    
    def process_metadata(self, h5_path):
        x = a4.meta.meta(h5_path, ".*", meta_fields=["series_id"])
        x["samples"] = x.index
        x = x[['series_id', 'samples']].drop_duplicates()
        return x

    def get_top_samples(self, df, n = 1000):
        df_meta = self.transcriptome_embedding.embeddings_index.copy()
        samps = df.iloc[:,1].sort_values(ascending = False).head(n).index.to_list()
        series_of_interest = df_meta[df_meta.index.isin(samps)]
        return series_of_interest

    def transcriptome_search(self, geneset, nsamples = 1000):
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
        series_of_interest = self.get_top_samples(final_df1, n = nsamples)
        relevant_series = series_of_interest["series_id"] # extract only studies (series)
        series_of_interest = series_of_interest.rename(columns={"series_id": "gse_id"})
        return  series_of_interest, relevant_series
        

    def semantic_search(self, test_text, k = 50):
        df_res = self.rag_embedding.query_rag(test_text, k = k)
        return df_res, df_res["gse_id"]

    def get_transcriptome_series_of_relevance_from_series(self, series_of_interest, k = 5):
        """use this in the series -> series on transcription
        """
        series_list = []
        target_series = series_of_interest
        for i in tqdm(target_series):
            res_temp = self.transcriptome_embedding.get_closest_transcriptional_studies(i, k = k)
            series_list.append(res_temp)
        additional_series = pd.concat(series_list)
        additional_series["gse_id_text"] = [self.rag_embedding.get_text_linked_to_gse(i) for i in additional_series["gse_id"]]
        return additional_series
    
    def get_semantic_series_of_relevance_from_series(self, series_of_interest, k = 5):
        """use this in the series -> series on semantic
        """
        additional_series = pd.concat([self.rag_embedding.get_closest_semantic_studies(i, k = k) for i in series_of_interest])
        return additional_series

    def transcriptome_search_with_semantic_expansion(self, geneset_query, search = 1000, expand= 5):
        if expand == 0:
            series_df, series_of_interest = self.transcriptome_search(geneset_query, search)
            return None, series_df
        else:
            series_df, series_of_interest = self.transcriptome_search(geneset=geneset_query, nsamples=search)
            additional_series = self.get_semantic_series_of_relevance_from_series(series_of_interest, expand)
            return additional_series, series_df

    def transcriptome_search_with_transcriptome_expansion(self, geneset_query, search = 1000, expand= 5):
        if expand == 0:
            series_df, series_of_interest = self.transcriptome_search(geneset_query, search)
            return None, series_df
        else:
            series_df, series_of_interest = self.transcriptome_search(geneset_query, search)
            additional_series = self.get_transcriptome_series_of_relevance_from_series(series_of_interest, expand)
            return additional_series, series_df

    def semantic_search_with_semantic_expansion(self, text_query, search = 50, expand= 5):
        print(text_query)
        if expand == 0:
            series_df, series_of_interest = self.semantic_search(text_query, search)
            return None, series_df
        else:
            series_df, series_of_interest = self.semantic_search(text_query, search)
            additional_series = self.get_semantic_series_of_relevance_from_series(series_of_interest, expand)
            return additional_series, series_df

    def semantic_search_with_transcriptome_expansion(self, text_query, search = 50, expand = 5):
        print(text_query)
        if expand == 0:
            series_df, series_of_interest = self.semantic_search(text_query, search)
            return None, series_df
        else:
            series_df, series_of_interest = self.semantic_search(text_query, search)
            additional_series = self.get_transcriptome_series_of_relevance_from_series(series_of_interest, expand)
            return additional_series, series_df

    def search(self, geneset, text_query, search = "semantic", expand = "transcriptome", perform_enrichment = False):
        """this is the main function
        """

        results_object = Results(None, None, None) # new results object

        if text_query is None:
            search = "transcriptome"

        if geneset is None:
            search = "semantic"

        if search == "semantic":
            if expand == "transcriptome":
                additional_series, seed_series = self.semantic_search_with_transcriptome_expansion(text_query)
                results_object.seed_studies = seed_series
                results_object.expansion_studies = additional_series
            elif expand == "semantic":
                additional_series, seed_series = self.semantic_search_with_semantic_expansion(text_query)
                results_object.seed_studies = seed_series
                results_object.expansion_studies = additional_series
        if search == "transcriptome":
            if expand == "transcriptome":
                additional_series, seed_series = self.transcriptome_search_with_transcriptome_expansion(geneset)
                results_object.seed_studies = seed_series
                results_object.expansion_studies = additional_series
            elif expand == "semantic":
                additional_series, seed_series = self.transcriptome_search_with_semantic_expansion(geneset)
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
    def __init__(self, seed_df, expansion_df, samples_df):
        self.seed_studies = seed_df
        self.expansion_studies = expansion_df
        self.samples = samples_df
    



        

