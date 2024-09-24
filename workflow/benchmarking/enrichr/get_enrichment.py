import os
import pandas as pd
from tqdm import tqdm
import traceback
import logging
from sample_explorer.sample_explorer import RNASeqAnalysis
from sample_explorer.utils import MsigDB_store
import io
import sys

# This script performs the following tasks:
# 1. Loads metadata and gene set information from various files
# 2. Initializes necessary objects for ssGSEA (using the RNAseqAnalysis module) and gene set handling
# 3. Processes each file in the results_enrichr_metadata directory
# 4. For each file, it:
#    a. Checks if the enrichment analysis has already been performed
#    b. Loads the corresponding gene set
#    c. Identifies relevant samples from the metadata
#    d. Performs enrichment analysis on the samples using the gene set
#    e. Saves the enrichment results to a CSV file
# 5. Logs the progress and any errors encountered during processing
# 6. Implements error handling and logging for robustness

# create "tempfiles/enrichr_anndata/"
os.makedirs("tempfiles/enrichr_anndata/", exist_ok=True)

def setup_logging(log_buffer=None):
    if log_buffer is None:
        log_buffer = io.StringIO()
    
    logging.basicConfig(
        stream=log_buffer,
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO
    )
    logging.info("Started processing.")
    return log_buffer

def load_data():
    res_files = os.listdir("tempfiles/results_enrichr_metadata/")
    x = RNASeqAnalysis("data/human_gene_v2.2.h5")
    df_series = pd.read_pickle("data/whole_metadata_human.p")
    msigdb_store = MsigDB_store("data/msigdb.v2023.2.Hs.symbols.gmt", "data/test_set_misgdb.csv")
    logging.info("Loaded data and initialized objects.")
    return res_files, x, df_series, msigdb_store

def process_file(filename, df_series, msigdb_store, x):
    dest_filename = "tempfiles/enrichr_anndata/" + filename.replace(".csv", "") + "_enrichment.csv"
    
    if os.path.exists(dest_filename):
        logging.info(f"File '{dest_filename}' exists, skipping.")
        return

    genesetname = filename.replace(".csv", "")
    geneset = msigdb_store.get_gene_set_by_name(genesetname)
    
    df = pd.read_csv(f"tempfiles/results_enrichr_metadata/{filename}", index_col=0)
    samps = df_series[df_series["series_id"].isin(df.index)].index.drop_duplicates()

    if len(samps) < 5000:
        anndata_res = x.perform_enrichment_on_samples(samps, geneset["geneset"])
    else:
        anndata_res = x.perform_enrichment_on_samples_batched(samps, geneset["geneset"])

    anndata_res.to_csv(dest_filename)
    logging.info(f"Successfully processed {filename} and saved to {dest_filename}.")

def main():
    log_buffer = setup_logging()
    
    res_files, x, df_series, msigdb_store = load_data()

    for filename in tqdm(res_files):
        try:
            process_file(filename, df_series, msigdb_store, x)
        except Exception as e:
            genesetname = filename.replace(".csv", "")
            logging.error(f"Failure processing {genesetname}: {str(e)}")
            logging.error(traceback.format_exc())
    
    with open("logfiles/processing_log.txt", "w") as log_file:
        log_file.write(log_buffer.getvalue())
    
    logging.info("Processing completed.")

if __name__ == "__main__":
    main()