import os
import pandas as pd
import pandas as pd
from tqdm import tqdm
import time
import os
from GEOparse import GEOparse
import signal

# Specify the folder path containing the Pickled files
folder_path = 'superseries'

# List all files in the folder with a .pkl extension
pickle_files = [file for file in os.listdir(folder_path) if file.endswith('.p')]

# Initialize an empty DataFrame to store the merged data
merged_df = []

# Loop through each Pickled file, read it, and append to the merged DataFrame
for file in pickle_files:
    file_path = os.path.join(folder_path, file)
    #print(file_path)
    df = pd.read_pickle(file_path).T
    merged_df.append(df)

def handler(signum, frame):
    raise TimeoutError("Time limit exceeded")
    
# Print the merged DataFrame or perform further operations

merged = pd.concat(merged_df)

currently_downloaded = merged.index.unique()

series_id = pd.read_csv("results/superseries_id_only.csv")["value"]

list_of_accessions = series_id[~series_id.isin(currently_downloaded)] # 84 studies not represented due to download errors

merged.to_csv("results/superseries.csv")