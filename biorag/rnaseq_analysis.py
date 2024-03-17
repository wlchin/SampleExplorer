
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


class RNASeqAnalysis:
    def __init__(self, file):
        self.file = file

    def create_anndata_from_series(self, series):
        series_counts = a4.data.series(self.file, series)
        metadata = a4.meta.series(self.file, series)
        test_ad = ad.AnnData(series_counts.T, dtype = np.float32)
        test_ad.obs = metadata
        test_ad.var_names_make_unique()
        return test_ad

    def create_anndata_from_samples(self, samples):
        sample_counts = a4.data.samples(self.file, samples)
        metadata = a4.meta.samples(self.file, samples)
        test_ad = ad.AnnData(sample_counts.T, dtype = np.float32)
        test_ad.obs = metadata
        test_ad.var_names_make_unique()
        return test_ad
    
    def list_to_dc_geneset(self, list_test_geneset):
        df = pd.DataFrame(list_test_geneset)
        df.columns = ["genesymbol"]
        df["geneset"] = "query_gene_set"
        return df

    def create_anndata_series(self, query):
        anndata_list = []

        for i in query:
            try:
                current_ad = self.create_anndata_from_series(i)
                anndata_list.append(current_ad)

            except Exception as e:
                
                print(f"An error occurred for element {i}: {e}")

        concat_adata = ad.concat(anndata_list, axis=0)
        return concat_adata

    def perform_enrichment_on_series(self, query_series, gene_set):
        
        current_ad = self.create_anndata_series(query_series)
        
        gene_set = self.list_to_dc_geneset(gene_set)
        
        dc.run_ora(
            mat=current_ad,
            net=gene_set,
            source='geneset',
            target='genesymbol',
            verbose=True, use_raw=False
        )

        results = current_ad.obsm['ora_estimate']
        results["pvals"] = current_ad.obsm['ora_pvals']

        df_res = pd.concat([results, current_ad.obs], axis=1)
        #max_index = df_res.groupby('series_id')['query_gene_set'].idxmax()

        #result_df = df_res.loc[max_index]
        return df_res

    def perform_enrichment_on_samples(self, query_sample_list, gene_set):
        
        current_ad = self.create_anndata_from_samples(query_sample_list)
        
        gene_set = self.list_to_dc_geneset(gene_set)
        
        dc.run_ora(
            mat=current_ad,
            net=gene_set,
            source='geneset',
            target='genesymbol',
            verbose=True, use_raw=False
        )

        results = current_ad.obsm['ora_estimate']
        results["pvals"] = current_ad.obsm['ora_pvals']

        df_res = pd.concat([results, current_ad.obs], axis=1)

        return df_res


    def perform_enrichment_on_samples_batched(self, element_list, geneset, batch_size = 5000):

        # Calculate the number of batches
        num_batches = (len(element_list) + batch_size - 1) // batch_size
        
        res_list = []
        
        # Iterate over batches
        for i in range(num_batches):
            start_index = i * batch_size
            end_index = (i + 1) * batch_size
            batch = element_list[start_index:end_index]
        
            # Apply your function to the current batch
            result = self.perform_enrichment_on_samples(batch, geneset)
            res_list.append(result)
        
        comb_df = pd.concat(res_list)
        return comb_df