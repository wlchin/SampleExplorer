import archs4py as a4
import pandas as pd

file = "human_gene_v2.2.h5"
df = a4.data.rand(file, 100, remove_sc=True)

x = pd.DataFrame(df.index)
x.index = x[0]
x.index.name = "gene_name"
x.columns = ["gene_name"]

x.to_pickle("results/genes_v2.p")