import os
import traceback
import pandas as pd
from tqdm import tqdm
import pickle
from biorag.biorag import (
    Transcriptome_embedding
)
from biorag.utils import MsigDB_store

# Purpose of script
# This script analyzes enrichment results for gene sets and calculates interpoint distances.
# It processes gene sets from MsigDB, loads or computes enrichment results,
# and calculates transcriptome-based distances between samples.
# Results are saved for each gene set and collected for further analysis.


# Initialize objects
import pickle

# Load pickle files
with open("data/important_data_for_LLM/transcription_index_v2.p", "rb") as f:
    transcription_index = pickle.load(f)

with open("data/important_data_for_LLM/transcription_embedding_matrix_v2.pkl", "rb") as f:
    transcription_embedding_matrix = pickle.load(f)

# Initialize objects
exp = Transcriptome_embedding(transcription_index, transcription_embedding_matrix)


msigdb_store = MsigDB_store(
    "data/msigdb.v2023.2.Hs.symbols.gmt", 
    "data/test_set_misgdb.csv"
)

gene_lists_collection = msigdb_store.list_genesets()

output_dir = "tempfiles/results_enrichr_distance/"
os.makedirs(output_dir, exist_ok=True)

# Processing loop
genesets = os.listdir("tempfiles/results_enrichr_metadata/") 

# remove the extension name from the strings in the list
filtered_list = [i.replace(".csv", "") for i in genesets]

results  = []

for i in tqdm(filtered_list):
    try:
        geneset_name = i
        output_filename = f"tempfiles/results_enrichr_distance/{geneset_name}_distance.pkl"
        
        # Check if the output file already exists
        if os.path.exists(output_filename):
            print(f"Loading existing results for {geneset_name}.")
            with open(output_filename, 'rb') as f:
                distance = pickle.load(f)
        else:
            # Load enrichment results from CSV
            enrichment_filename = f"tempfiles/enrichr_anndata/{geneset_name}_enrichment.csv"
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
results_df.to_csv("tempfiles/enrichr_collected_distances.csv", index=False)

print("Processing complete, results saved to 'results/enrichr_collected_distances.csv'")
