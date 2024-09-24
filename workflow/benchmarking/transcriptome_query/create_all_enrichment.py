from tqdm import tqdm
import pandas as pd
from biorag.utils import MsigDB_store 
from biorag.biorag import RNASeqAnalysis, Transcriptome_enrichment
from biorag.enrichment import Transcriptome_enrichment
import traceback
from typing import List
import anndata as ad
import os

# creates the enrichment df

#The transcriptome enrichment modiule needs to accept the anndata object. 
#Hence you have to load it, and then pass it to the function.

import anndata as ad

condition_matrix = ad.read("data/test_transcriptome_db.h5ad")

te = Transcriptome_enrichment(condition_matrix)

x = RNASeqAnalysis("data/human_gene_v2.2.h5")

msigdb_store = MsigDB_store("data/msigdb.v2023.2.Hs.symbols.gmt", "data/test_set_misgdb.csv")

DATA_DIR = "data/"
GENESET_DIR = "gene_sets/"
TARGET_DIR = "tempfiles/enrichment/"

os.makedirs(TARGET_DIR, exist_ok=True)

def load_test_genesets() -> List[str]:
    test_gs = pd.read_csv(f"{GENESET_DIR}test_gene_sets.txt", header=None)[0].to_list()
    return test_gs

my_list = te.memmap_adata.obs.index # create the memory map for the conditions - there are 287,000 conditions

testset = load_test_genesets()

dict_genes = {}

for i in testset:
    res = msigdb_store.get_gene_set_by_name(i)
    dict_genes[i] = res["geneset"]

user_batch_size = 250

output_df1_list = []
output_df2_list = []

# Iterate over the list in batches with tqdm
for i in tqdm(range(0, len(my_list), user_batch_size), desc="Processing Batches", unit="batch"):
    try:
        samps = my_list[i:i + user_batch_size]
        df1, df2 = te.run_decouplr_on_memmaped_adata_with_samples(dict_genes, samps)

        # Store the output from each iteration
        output_df1_list.append(df1)
        output_df2_list.append(df2)
        
        if i % 1000 == 0:
            temp_df1 = pd.concat(output_df1_list).astype('float32')
            temp_df2 = pd.concat(output_df2_list).astype('float32')
            
            fileval = "tempfiles/enrichment/pval_df_" + str(i) + ".p"
            fileval2 = "tempfiles/enrichment/enrichment_df_" + str(i) + ".p"

            temp_df1.to_pickle(fileval)
            temp_df2.to_pickle(fileval2)
            
            output_df1_list_temp = []
            output_df2_list_temp = []
            
    except:
        print("failure", "batch: ", str(i))
        print(traceback.print_exc())
    
final_df1 = pd.concat(output_df1_list).astype('float32')
final_df2 = pd.concat(output_df2_list).astype('float32')

final_df1.to_pickle("tempfiles/enrichment/final_pval_df.p")
final_df2.to_pickle("tempfiles/enrichment/final_enrichment_df.p")