import pandas as pd
import os
import pandas as pd

x = pd.read_csv("data/gene_sets_msigdb.csv")
# remove those rows with NaN values for description_full and description_brief
x = x.dropna(subset=["description_full", "description_brief"])
x.to_csv("data/test_set_misgdb.csv", index=False)
x.sample(1000)["standard_name"].to_csv("data/test_gene_sets.txt", index=False, header=False)




