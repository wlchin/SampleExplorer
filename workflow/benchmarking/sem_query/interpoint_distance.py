import os
import traceback
import pandas as pd
from tqdm import tqdm
import pickle
from sample_explorer.sample_explorer import (
    Transcriptome_embedding
)
from sample_explorer.utils import MsigDB_store
from typing import List

# Purpose of script
# This script analyzes enrichment results for gene sets and calculates interpoint distances.
# It processes gene sets from MsigDB, loads or computes enrichment results,
# and calculates transcriptome-based distances between samples.
# Results are saved for each gene set and collected for further analysis.

# Initialize objects
import pickle


DATA_DIR = "data/"
QUERY_DIR = "tempfiles/sem_and_tran_df/"
TARGET_DIR = "tempfiles/sem_and_tran_df_distance/"
GENESET_DIR = "gene_sets/"

def load_test_genesets() -> List[str]:
    test_gs = pd.read_csv(f"{GENESET_DIR}test_gene_sets.txt", header=None)[0].to_list()
    return test_gs

# Load pickle files
def create_transcriptome_embedding():
    with open("data/transcription_index_v2.p", "rb") as f:
        transcription_index = pickle.load(f)

    with open("data/transcription_embedding_matrix_v2.pkl", "rb") as f:
        transcription_embedding_matrix = pickle.load(f)

# Initialize objects
    exp = Transcriptome_embedding(transcription_index, transcription_embedding_matrix)
    return exp


def check_folder_existence():
    os.makedirs(TARGET_DIR, exist_ok=True)


exp = create_transcriptome_embedding()
msigdb_store = MsigDB_store(
    "data/msigdb.v2023.2.Hs.symbols.gmt", 
    "data/test_set_misgdb.csv"
)
gene_lists_collection = msigdb_store.list_genesets()
# Processing loop

filtered_list = load_test_genesets()

results  = []
check_folder_existence()

for i in tqdm(filtered_list):
    try:
        geneset_name = i
        output_filename = f"{TARGET_DIR}{geneset_name}_distance.pkl"
        
        # Check if the output file already exists
        if os.path.exists(output_filename):
            print(f"Loading existing results for {geneset_name}.")
            with open(output_filename, 'rb') as f:
                distance = pickle.load(f)
        else:
            # Load enrichment results from CSV
            enrichment_filename = f"{QUERY_DIR}{geneset_name}.csv"
            enrichment_df = pd.read_csv(enrichment_filename, index_col=0)
            
            # Extract samples from the index of the dataframe
            samples = enrichment_df.index.tolist()
            
            # Calculate interpoint distance for the samples
            distance = exp.calculate_interpoint_distance_for_samples(samples, p=1)
            
            # Save the result to a file
            with open(output_filename, 'wb') as f:
                pickle.dump(distance, f)
        
        # Collect results into a list of dictionaries
        results.append({
            "geneset_name": geneset_name,
            "distance": distance
        })

    except Exception as e:
        print(f"Failure: {geneset_name}")
        print(traceback.format_exc())

# Convert the results list into a DataFrame
results_df = pd.DataFrame(results)

# Save the DataFrame to a CSV file
results_df.to_csv("tempfiles/sem_and_tran_collected_distances.csv", index=False)

print("Processing complete, results saved to tempfiles/sem_and_tran_collected_distances.csv")
