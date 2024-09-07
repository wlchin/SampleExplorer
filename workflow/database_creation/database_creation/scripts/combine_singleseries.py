import os
import pandas as pd

# Specify the folder path containing the Pickled files
folder_path = 'single_series'

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

# Print the merged DataFrame or perform further operations

merged = pd.concat(merged_df)

merged.to_csv("results/single_series.csv")