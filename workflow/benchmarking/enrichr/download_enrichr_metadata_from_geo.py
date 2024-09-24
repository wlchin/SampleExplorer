import os
import re
import time
import traceback
import pandas as pd
import GEOparse
from tqdm import tqdm
from io import StringIO

# This script downloads and processes metadata from GEO (Gene Expression Omnibus) for Enrichr analysis.
# It performs the following tasks:
# 1. Reads a list of GEO accession IDs from input files
# 2. Downloads metadata for each accession using GEOparse
# 3. Collects and organizes the metadata into a DataFrame
# 4. Handles errors and logs issues during the download process
# 5. Saves the collected metadata to a CSV file for further analysis
# 6. Implements rate limiting to avoid overloading the GEO server


def collect_studies(list_of_accessions, log_buffer):
    """
    Collect metadata for a list of GEO accessions.

    Args:
        list_of_accessions (list): List of GEO accession IDs.
        log_buffer (StringIO): Buffer to store log messages.

    Returns:
        pd.DataFrame: DataFrame containing the metadata of the GEO accessions.
    """
    dict_to_df = {}
    
    for accession in tqdm(list_of_accessions, desc="Processing GEO Accessions"):
        try:
            gse = GEOparse.get_GEO(geo=accession, destdir="./", include_data=False, silent=True)
            time.sleep(1)
            file_to_remove = f"{accession}_family.soft.gz"
            dict_to_df[accession] = gse.metadata
    
            # Check and remove the downloaded file if it exists
            if os.path.exists(file_to_remove):
                os.remove(file_to_remove)
        
        except Exception as e:
            error_message = f"Error processing accession {accession}: {e}"
            print(error_message)
            log_buffer.write(error_message + "\n")
            log_buffer.write(traceback.format_exc() + "\n")
    
    geo_accession_df = pd.DataFrame(dict_to_df).T
    return geo_accession_df

def grab_study_GEO(input_string):
    """
    Extract the first GEO accession (GSEXXX) from a string.

    Args:
        input_string (str): Input string containing the GEO accession.

    Returns:
        str: The GEO accession found in the input string, or None if not found.
    """
    pattern = r'GSE\d+'
    match = re.search(pattern, input_string)
    
    if match:
        return match.group()
    else:
        print(f"No GEO accession found in: {input_string}")
        return None

def process_enrichr_files(enrichr_dir, metadata_dir, log_buffer, p_value_threshold=0.05):
    """
    Process Enrichr results files, extract significant studies, and collect metadata.

    Args:
        enrichr_dir (str): Directory containing Enrichr result files.
        metadata_dir (str): Directory to save metadata results.
        log_buffer (StringIO): Buffer to store log messages.
        p_value_threshold (float): Threshold for selecting significant studies.
    """
    os.makedirs(metadata_dir, exist_ok=True)
    enrichr_files = os.listdir(enrichr_dir)

    for filename in tqdm(enrichr_files, desc="Processing Enrichr Files"):
        try:
            load_file = os.path.join(enrichr_dir, filename)
            savefile_metadata = os.path.join(metadata_dir, filename)
            
            # Load Enrichr results
            enrichr_results = pd.read_csv(load_file)
            
            # Filter significant studies
            sig_studies = enrichr_results[enrichr_results["Adjusted P-value"] < p_value_threshold]["Term"]
            sig_study_list = [grab_study_GEO(term) for term in sig_studies if grab_study_GEO(term) is not None]
            
            # Collect metadata for significant studies
            metadata_df = collect_studies(sig_study_list, log_buffer)
            metadata_df.to_csv(savefile_metadata)
            
            log_buffer.write(f"Successfully processed file: {filename}\n")
        
        except Exception as e:
            error_message = f"Error processing file {filename}: {e}"
            print(error_message)
            log_buffer.write(error_message + "\n")
            log_buffer.write(traceback.format_exc() + "\n")


log_buffer = StringIO()

enrichr_dir = "tempfiles/results_enrichr"
metadata_dir = "tempfiles/results_enrichr_metadata"

process_enrichr_files(enrichr_dir, metadata_dir, log_buffer)

# Write the accumulated log information to a log file at the end
with open("logfiles/get_metadata_enrichr.log", "w") as log_file:
    log_file.write(log_buffer.getvalue())