
import pandas as pd
import decoupler as dc
import numpy as np
import statsmodels.stats.multitest as smm

class Transcriptome_enrichment:
    def __init__(self, h5ad_obj):
        """
        Initializes the Enrichment class, which stores representative transcriptomes

        Parameters:
        - h5ad_obj: An h5ad object representing the data to be enriched.
        """
        self.memmap_adata = h5ad_obj
   
    def list_to_dc_geneset(self, list_test_geneset):
        """
        Convert a list of genes to a pandas DataFrame representing a gene set.

        Args:
            list_test_geneset (list): A list of gene symbols.

        Returns:
            pandas.DataFrame: A DataFrame with two columns: 'genesymbol' and 'geneset'.
                              The 'genesymbol' column contains the gene symbols from the input list,
                              and the 'geneset' column is set to 'query_gene_set' for all rows.

        """
        df = pd.DataFrame(list_test_geneset)
        df.columns = ["genesymbol"]
        df["geneset"] = "query_gene_set"
        return df
    
    def list_to_dc_geneset_dictionary(self, list_test_geneset_dict):
        """
        Convert a list of test geneset dictionaries into a pandas DataFrame.

        Args:
            list_test_geneset_dict (list): A list of test geneset dictionaries.

        Returns:
            pandas.DataFrame: A DataFrame containing the genesets and their corresponding genes.

        """
        geneset_df_list = []
        for i, j in list_test_geneset_dict.items():
            df = pd.DataFrame(j)
            df.columns = ["genesymbol"]
            df["geneset"] = i
            geneset_df_list.append(df)
        return pd.concat(geneset_df_list)
    
    def run_decouplr_on_memmaped_adata_with_samples(self, genelist, samples_list):
        """
        Runs the decouplr algorithm on a memmaped AnnData object with a specified list of samples.

        Args:
            genelist (dict or list): The gene list to use for the decouplr algorithm. It can be either a dictionary
                representing a geneset or a list of genes.
            samples_list (list): The list of sample names to include in the analysis.

        Returns:
            pandas.DataFrame: A DataFrame containing the results of the decouplr algorithm, including the p-values,
            estimates, and corrected p-values (if applicable).
        """
        if isinstance(genelist, dict):
            geneset = self.list_to_dc_geneset_dictionary(genelist)
        elif isinstance(genelist, list):
            geneset = self.list_to_dc_geneset(genelist)
        else:
            print(genelist)

        adata = self.memmap_adata
        query_adata = adata[adata.obs.index.isin(samples_list)].to_memory()
        query_adata.var_names_make_unique()
        np.nan_to_num(query_adata.X, copy=False, nan=0.0, posinf=0.0, neginf=0.0)

        dc.run_ora(
            mat=query_adata,
            net=geneset,
            source='geneset',
            target='genesymbol',
            verbose=True,
            use_raw=False,
            min_n=1
        )

        if query_adata.obsm["ora_pvals"].shape[1] > 1:
            return query_adata.obsm["ora_pvals"], query_adata.obsm["ora_estimate"]
        else:
            res = pd.concat([query_adata.obsm["ora_pvals"], query_adata.obsm["ora_estimate"]], axis=1)
            res.columns = ["pvals", "estimate"]
            _, p_corrected, _, _ = smm.multipletests(res['pvals'], method='fdr_bh')
            res['p_corrected'] = p_corrected
            return res
    

