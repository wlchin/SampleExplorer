import pandas as pd
from tqdm import tqdm
import time
import os
from GEOparse import GEOparse

save_interval = 100  # Set the save interval

df_melted = pd.read_csv("results/no_superseries_ids_only.csv")

list_of_accessions = df_melted["gse_id"]
dict_to_df = {}

save_folder = "single_series"
# Create the folder if it doesn't exist
if not os.path.exists(save_folder):
    os.makedirs(save_folder)
# Set the save interval  # Replace "destination_folder" with the desired folder path

dict_to_df = {}  # Initialize an empty dictionary

for idx, i in enumerate(tqdm(list_of_accessions)):
    try:
        gse = GEOparse.get_GEO(geo=i, destdir="./", include_data=False, silent=True)
        time.sleep(3)
        file_to_remove = i + "_family.soft.gz"
        dict_to_df[i] = gse.metadata

        # Check if the file exists before attempting to remove it
        if os.path.exists(file_to_remove):
            # Remove the file
            os.remove(file_to_remove)
        else:
            pass
        
        # Save DataFrame to disk every save_interval iterations
        if (idx + 1) % save_interval == 0:
            partial_df = pd.DataFrame(dict_to_df)
            partial_df.to_pickle(os.path.join(save_folder, f"metadata_{idx + 1}_new2.p"))
            
            # Clear the dictionary to free up memory
            dict_to_df = {}

    except Exception as e:
        print(f"Error processing accession {i}: {e}")
        pass

# Save the final DataFrame if there are remaining entries in the dictionary
if dict_to_df:
    final_df = pd.DataFrame(dict_to_df)
    final_df.to_pickle(os.path.join(save_folder, "metadata_final.p"))