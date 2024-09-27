
import pandas as pd
import anndata as ad
import decoupler as dc
import numpy as np
import archs4py as a4
import statsmodels.stats.multitest as smt
from tqdm import tqdm

class RNASeqAnalysis:
    def __init__(self, file):
        """
        Initializes a new instance of the `rnaseq_analysis` class.

        Args:
            file (str): The file path to be assigned to the `file` attribute.

        Returns:
            None
        """
        self.file = file

    def create_anndata_from_series(self, series):
        """Create an AnnData object with counts and metadata for a given series.

        Args:
            series (str): The name of the series.

        Returns:
            AnnData: An AnnData object containing the counts and metadata for the series.
        """
        series_counts = a4.data.series(self.file, series)
        metadata = a4.meta.series(self.file, series)
        test_ad = ad.AnnData(series_counts.T)
        test_ad.obs = metadata
        test_ad.var_names_make_unique()
        return test_ad

    def create_anndata_from_samples(self, samples):
        """Create an AnnData object with counts and metadata for the given samples.

        Args:
            samples (list): A list of sample names.

        Returns:
            AnnData: An AnnData object containing the counts and metadata for the samples.
        """
        sample_counts = a4.data.samples(self.file, samples)
        metadata = a4.meta.samples(self.file, samples)
        test_ad = ad.AnnData(sample_counts.T)
        test_ad.obs = metadata
        test_ad.var_names_make_unique()
        return test_ad
    
    def list_to_dc_geneset(self, list_test_geneset):
        """
        Convert a list of genes into a pandas DataFrame representing a geneset.

        Args:
            list_test_geneset (list): A list of genes.

        Returns:
            pandas.DataFrame: A DataFrame representing the geneset, with columns "genesymbol" and "geneset".

        """
        df = pd.DataFrame(list_test_geneset)
        df.columns = ["genesymbol"]
        df["geneset"] = "enrichment_score"
        return df

    def perform_enrichment_on_samples_batched(self, element_list, geneset, batch_size=2500):
        """
        Perform enrichment analysis on a list of elements in batches.

        Args:
            element_list (list): A list of elements to perform enrichment analysis on.
            geneset (str): The geneset to use for enrichment analysis.
            batch_size (int, optional): The size of each batch. Defaults to 5000.

        Returns:
            pandas.DataFrame: A combined DataFrame containing the results of enrichment analysis for all batches.
        """

        # Calculate the number of batches
        num_batches = (len(element_list) + batch_size - 1) // batch_size

        res_list = []

        # Iterate over batches
        for i in tqdm(range(num_batches)):
            start_index = i * batch_size
            end_index = (i + 1) * batch_size
            batch = element_list[start_index:end_index]

            # Apply your function to the current batch
            result = self.perform_enrichment_on_samples(batch, geneset)
            res_list.append(result)

        comb_df = pd.concat(res_list)
        return comb_df

    def perform_enrichment_on_series(self, query_series, gene_set):
        """
        This function accepts a list of series and a gene set and returns a DataFrame with the enrichment results.

        Parameters:
            query_series (list): A list of series.
            gene_set (str): The gene set to perform enrichment on.

        Returns:
            pandas.DataFrame: A DataFrame with the enrichment results.
        """
        
        samps = self.create_samples_from_series(query_series)
        samples_metadata = pd.concat(samps) # concatenate dataframes
        df_res = self.perform_enrichment_on_samples_batched(samples_metadata.index.drop_duplicates(), gene_set)
        return df_res
    
    def create_samples_from_series(self, series):
        """
        Create samples from a given series.

        Parameters:
        - series: The series to create samples from.

        Returns:
        - (list of ) sample dataframes: The samples created from the series.
        """
        samples = [a4.meta.series(self.file, i) for i in series]
        return samples

    def perform_enrichment_on_samples(self, query_sample_list, gene_set): 
        """
        Perform enrichment analysis on a list of samples using a given gene set.

        Args:
            query_sample_list (list): A list of samples to perform enrichment analysis on.
            gene_set (list): A list of genes representing the gene set to use for enrichment analysis.

        Returns:
            pandas.DataFrame: A DataFrame containing the enrichment results.

        """
        
        current_ad = self.create_anndata_from_samples(query_sample_list)
        
        gene_set = self.list_to_dc_geneset(gene_set)
        
        dc.run_ora(
            mat=current_ad,
            net=gene_set,
            source='geneset',
            target='genesymbol',
            verbose=False, use_raw=False, min_n = 1,
        )

        results = current_ad.obsm['ora_estimate']
        results["pvals"] = current_ad.obsm['ora_pvals']
        results["p_corrected"] = smt.multipletests(results["pvals"], method='fdr_bh')[1]

        df_res = pd.concat([results, current_ad.obs], axis=1)

        return df_res


