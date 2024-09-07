import archs4py as a4
import pandas as pd

file = "human_gene_v2.2.h5"

metadata = a4.meta.meta(file, ".*", meta_fields=["series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"])
metdata_clean = metadata.drop_duplicates()
metdata_clean.to_csv("results/metadata_geo_human.csv")

result_superseries = metdata_clean[metdata_clean['series_id'].str.contains(',')]
result_superseries.to_csv("results/arch_geo_superseries.csv")

result_no_superseries = metdata_clean[~(metdata_clean['series_id'].str.contains(','))]
result_no_superseries.to_csv("results/arch_geo_nosuperseries.csv")

x = pd.read_csv("results/arch_geo_nosuperseries.csv")
x1 = pd.DataFrame(x["series_id"].unique())
x1.columns = ["gse_id"]
x1.to_csv("results/no_superseries_ids_only.csv")

x2 = pd.read_csv("results/arch_geo_superseries.csv")
x2_parsed = pd.DataFrame(x2["series_id"].str.split(",", expand = True).stack().reset_index(drop=True).unique()).dropna()
x2_parsed.columns = ["gse_id"]# 3106
x2_parsed.to_csv("results/superseries_ids_only.csv")