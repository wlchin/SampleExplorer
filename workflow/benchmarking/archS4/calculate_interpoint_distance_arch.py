
import os
import traceback
import pandas as pd
import pickle
from sample_explorer.sample_explorer import (
    Transcriptome_embedding
)
from sample_explorer.utils import MsigDB_store

# Purpose of script
# This script analyzes ARCHS4 results for gene sets and calculates interpoint distances.
# It processes gene sets from MsigDB, loads ARCHS4 results from pickle files,
# and calculates transcriptome-based distances between samples.
# Results are saved for each gene set and collected for further analysis.

# Initialize objects
import pickle

# Load pickle files
with open("data/transcription_index_v2.p", "rb") as f:
    transcription_index = pickle.load(f)

with open("data/transcription_embedding_matrix_v2.pkl", "rb") as f:
    transcription_embedding_matrix = pickle.load(f)

# Initialize objects
exp = Transcriptome_embedding(transcription_index, transcription_embedding_matrix)


msigdb_store = MsigDB_store(
    "data/msigdb.v2023.2.Hs.symbols.gmt", 
    "data/test_set_misgdb.csv"
)

BASE_DIR = "tempfiles/arch_enrichment/"
OUTPUT_DIR = "tempfiles/results_archs4_distance/"

os.makedirs(OUTPUT_DIR, exist_ok=True)

files = os.listdir(BASE_DIR)

results = []

for i in files:
    try:
        geneset_name = i.replace(".csv", "")
        output_filename = f"{OUTPUT_DIR}{geneset_name}_distance.pkl"
        
        # Check if the output file already exists
        if os.path.exists(output_filename):
            print(f"Loading existing results for {geneset_name}.")
            with open(output_filename, 'rb') as f:
                distance = pickle.load(f)
        else:
            # Load ARCHS4 results from pickle file
            archs4_filename = f"{BASE_DIR}/{geneset_name}.csv"
            archs4_df = pd.read_csv(archs4_filename, index_col=0)   
            
            # Extract top 1000 samples based on Z-score
            samples = archs4_df.index # implement filtering if required
            
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
results_df.to_csv("tempfiles/archs4_collected_distances.csv", index=False)

print("Processing complete, results saved to tempfiles/archs4_collected_distances.csv")


