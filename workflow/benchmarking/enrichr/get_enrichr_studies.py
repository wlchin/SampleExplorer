import os
import time
import traceback
import gseapy as gp
from io import StringIO
from biorag.utils import MsigDB_store
import pandas as pd


# This script performs the following tasks:
# 1. Loads gene sets from MsigDB and filters for those containing "UP"
# 2. Retrieves Enrichr results for each filtered gene set using the RNAseq_Automatic_GEO_Signatures_Human_Up database
# 3. Saves the Enrichr results for each gene set as a pickle file
# 4. Implements error handling and logging for robustness
# 5. Uses rate limiting to avoid overloading the Enrichr server


# Initialize the MsigDB store
msigdb_store = MsigDB_store("data/msigdb.v2023.2.Hs.symbols.gmt", "data/test_set_misgdb.csv")
gene_lists_collection = msigdb_store.list_genesets()

# Filter gene sets containing "UP"
#filtered_gene_sets = [gene_set for gene_set in gene_lists_collection if "UP" in gene_set]
filtered_gene_sets = pd.read_csv("enrichr/test_data/test_gene_sets.txt", header=None)[0].to_list() # test set

# Define output directory and ensure it exists
output_dir = "tempfiles/results_enrichr/"
os.makedirs(output_dir, exist_ok=True)

# Initialize a StringIO buffer to accumulate log messages
log_buffer = StringIO()

def get_enrichr_result(gene_list):
    """
    Function to retrieve Enrichr results for a given gene list.
    """
    enrichment_result = gp.enrichr(
        gene_list=gene_list,
        gene_sets=['RNAseq_Automatic_GEO_Signatures_Human_Up'],
        organism='human',  # Set organism to human
        outdir=None  # Do not write results to disk
    )
    return enrichment_result.res2d

# Process each gene set
for gene_set_name in filtered_gene_sets:
    try:
        # Introduce a delay to avoid overloading the server
        time.sleep(1)

        # Retrieve the gene set from the MsigDB store
        gene_set = msigdb_store.get_gene_set_by_name(gene_set_name)["geneset"]

        # Get Enrichr results for the gene set
        enrichr_results = get_enrichr_result(gene_set)

        # Define the output filename
        output_file = os.path.join(output_dir, f"{gene_set_name}.csv")

        # Save the results to a pickle file
        enrichr_results.to_csv(output_file)

        # Log success
        log_buffer.write(f"Successfully processed gene set: {gene_set_name}\n")
    
    except Exception as e:
        error_message = f"Failed to process gene set '{gene_set_name}': {e}"
        log_buffer.write(error_message + "\n")
        
        # Capture and log the traceback
        log_buffer.write(traceback.format_exc() + "\n")
        
        # Save the traceback to an individual error file
        error_log_file = os.path.join(output_dir, f"{gene_set_name}_error.log")
        with open(error_log_file, "w") as ef:
            ef.write(traceback.format_exc())

# Write the accumulated log information to a log file at the end
with open("logfiles/enrichr_api.log", "w") as log_file:
    log_file.write(log_buffer.getvalue())