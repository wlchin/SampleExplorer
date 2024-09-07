import archs4py as a4
import numpy as np
from tqdm import tqdm

file = "human_gene_v2.2.h5"

meta_meta = a4.meta.meta(file, ".*", meta_fields=["series_id", "characteristics_ch1"])
df_no_duplicates = meta_meta.drop_duplicates()
df_no_duplicates.to_pickle("results/transcription_index_v2.pkl")
list_of_tuples = df_no_duplicates.apply(tuple, axis=1).tolist()

def get_average_counts(series_id, characteristics):
    samples = meta_meta[(meta_meta["series_id"] == series_id) & \
                    (meta_meta["characteristics_ch1"] == characteristics)].index.to_list()
    sample_counts = a4.data.samples(file, samples)
    av_counts = sample_counts.mean(1).values
    return av_counts

rows = 287553
cols = 67186

# Create a memmap array
memmap_filename = "data/study_counts.dat"
memmap_matrix = np.memmap(memmap_filename, dtype='float16', mode='w+', shape=(rows, cols))

# Number of rows to add
num_rows_to_add = 287553

# Sequentially add random rows
for i in tqdm(range(num_rows_to_add), desc="Adding Rows", unit="row"):
    study_id, characteristics = list_of_tuples[i]
    new_row = get_average_counts(study_id, characteristics)
    memmap_matrix[i, :] = new_row    
    memmap_matrix.flush()


memmap_matrix.flush()

# Close the memmap file
del memmap_matrix