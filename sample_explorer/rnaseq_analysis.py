
import pandas as pd
import anndata as ad
from anndata.utils import make_index_unique
import numpy as np
import logging
from sample_explorer import archs4
from sample_explorer.ora import bh_adjust, run_ora
from tqdm import tqdm

logger = logging.getLogger(__name__)


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
        return self.create_anndata_from_samples(archs4.samples_in_series(self.file, series))

    def create_anndata_from_samples(self, samples):
        """Create an AnnData object with counts and metadata for the given samples.

        Args:
            samples (list): A list of sample names.

        Returns:
            AnnData: An AnnData object containing the counts and metadata for the samples.
        """
        if isinstance(samples, str):
            raise TypeError("samples must be an iterable of sample IDs, not a string.")
        samples = list(samples)
        sample_counts = archs4.read_counts(self.file, samples)
        if sample_counts.shape[1] == 0:
            raise ValueError("None of the requested samples were found in the ARCHS4 file")
        missing = len(set(samples)) - sample_counts.shape[1]
        if missing:
            logger.warning("%d requested sample IDs not found in the ARCHS4 file", missing)
        metadata = archs4.read_sample_metadata(self.file, sample_counts.columns)
        # ARCHS4 repeats some gene symbols; dedupe as var_names_make_unique would, before AnnData warns
        sample_counts.index = make_index_unique(sample_counts.index)
        test_ad = ad.AnnData(sample_counts.T)
        test_ad.obs = metadata.loc[sample_counts.columns]
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
        samples = [archs4.read_sample_metadata(self.file, archs4.samples_in_series(self.file, i)) for i in series]
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
        
        pvals, scores = run_ora(current_ad.X, current_ad.obs_names, current_ad.var_names,
                                {"enrichment_score": list(gene_set)})
        results = pd.DataFrame({"enrichment_score": scores["enrichment_score"],
                                "pvals": pvals["enrichment_score"]})
        results["p_corrected"] = bh_adjust(results["pvals"])

        df_res = pd.concat([results, current_ad.obs], axis=1)

        return df_res


