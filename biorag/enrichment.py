
import pandas as pd
import decoupler as dc
import numpy as np
import statsmodels.stats.multitest as smm

class Transcriptome_enrichment:
    def __init__(self, h5ad_obj):
        self.memmap_adata = h5ad_obj
   
    def list_to_dc_geneset(self, list_test_geneset):
        df = pd.DataFrame(list_test_geneset)
        df.columns = ["genesymbol"]
        df["geneset"] = "query_gene_set"
        return df
    
    def list_to_dc_geneset_dictionary(self, list_test_geneset_dict):
        geneset_df_list = []
        for i,j in list_test_geneset_dict.items():
            df = pd.DataFrame(j)
            df.columns = ["genesymbol"]
            df["geneset"] = i
            geneset_df_list.append(df)
        return pd.concat(geneset_df_list)
    
    def run_decouplr_on_memmaped_adata_with_samples(self, genelist, samples_list):

        if isinstance(genelist, dict):
            geneset = self.list_to_dc_geneset_dictionary(genelist)
            #print(geneset)

        elif isinstance(genelist, list):
            geneset = self.list_to_dc_geneset(genelist)
            #print(geneset)

        else:
            print(genelist)

        adata = self.memmap_adata
        query_adata = adata[adata.obs.index.isin(samples_list)].to_memory()
        np.nan_to_num(query_adata.X, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        #query_adata.var_names_make_unique()

        dc.run_ora(
            mat=query_adata,
            net=geneset,
            source='geneset',
            target='genesymbol',
            verbose=True, use_raw=False
        )
        if query_adata.obsm["ora_pvals"].shape[1] > 1:
            return query_adata.obsm["ora_pvals"], query_adata.obsm["ora_estimate"]
        else:
            res = pd.concat([query_adata.obsm["ora_pvals"], query_adata.obsm["ora_estimate"]], axis = 1)
            res.columns = ["pvals", "estimate"]
            _, p_corrected, _, _ = smm.multipletests(res['pvals'], method='fdr_bh')
            res['p_corrected'] = p_corrected
            return res
    

