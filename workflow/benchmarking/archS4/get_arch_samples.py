import os
import pandas as pd
import logging
import traceback
from io import StringIO
from sample_explorer.sample_explorer import RNASeqAnalysis
from sample_explorer.utils import MsigDB_store, ARCHS4_API_query

#Initialize logging with a StringIO buffer
log_buffer = StringIO()

logging.basicConfig(
    stream=log_buffer,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# Initialize stores and analysis objects
data_source = "data/"
msigdb_store = MsigDB_store(data_source + "msigdb.v2023.2.Hs.symbols.gmt", data_source + "test_set_misgdb.csv")
archquery = ARCHS4_API_query(data_source + "whole_metadata_human.p")
rna_seq_analysis = RNASeqAnalysis(data_source + "human_gene_v2.2.h5")

# Filter gene sets containing "UP"
filtered_gene_sets = pd.read_csv("archs4/test_data/test_gene_sets.txt", header = None)[0].to_list()
#filtered_gene_sets = [gene_set for gene_set in msigdb_store.list_genesets() if "UP" in gene_set]

# Define output directory and ensure it exists
output_dir = "tempfiles/arch_enrichment/" # i will swtich this to arch_int_3
os.makedirs(output_dir, exist_ok=True)

# Process each gene set
for gene_set_name in filtered_gene_sets:
    try:
        # Generate the filename for the output
        output_file = os.path.join(output_dir, f"{gene_set_name}.csv")

        # Retrieve the gene set and associated genes
        gene_set = msigdb_store.get_gene_set_by_name(gene_set_name)["geneset"]

        # Query relevant series and samples
        _, samples = archquery.return_relevant_series_and_samples("up", gene_set)

        # Perform enrichment analysis on the samples
        enrichment_results = rna_seq_analysis.perform_enrichment_on_samples(samples, gene_set)

        # Save the results
        enrichment_results.to_csv(output_file)
        logging.info(f"Successfully processed gene set: {gene_set_name}")
    
    except Exception as e:
        error_message = f"Failed to process gene set '{gene_set_name}': {e}"
        logging.error(error_message)

        # Save traceback to a separate error file
        error_file = os.path.join(output_dir, f"{gene_set_name}_error.log")
        with open(error_file, "w") as ef:
            ef.write(traceback.format_exc())
        
        logging.info(f"Error details saved to {error_file}")

os.makedirs("logfiles/", exist_ok=True)
# Write the accumulated log information to a log file
with open("logfiles/arch_api.log", "w") as log_file:
    log_file.write(log_buffer.getvalue())