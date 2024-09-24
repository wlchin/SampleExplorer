import archs4py as a4
import numpy as np
from tqdm import tqdm
import pandas as pd

file = "data/human_gene_v2.2.h5"

meta_meta_no_sc = a4.meta.meta(file, ".*", meta_fields=["series_id", "characteristics_ch1"], remove_sc = True)
meta_meta_with_sc = a4.meta.meta(file, ".*", meta_fields=["series_id", "characteristics_ch1"])

df_no_duplicates_single_cell_excluded = meta_meta_no_sc.drop_duplicates()
df_no_duplicates_single_cell_excluded.to_pickle("data/no_single_cell.p")

df_no_duplicates_with_single_cell = meta_meta_with_sc.drop_duplicates()
df_no_duplicates_with_single_cell.to_pickle("data/whole_metadata_human.p")
