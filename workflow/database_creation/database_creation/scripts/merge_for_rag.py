import pandas as pd

comb = pd.read_csv("results/full_nosuperseries_with_pmid.csv", index_col=0)

comb.columns = ['Unnamed: 0', 'series_id', 'characteristics_ch1',
       'extract_protocol_ch1', 'source_name_ch1', 'title_ex', 'title', 'summary',
       'overall_design', 'pmid']

super1 = comb[['series_id', 'characteristics_ch1','title', 'summary',
       'overall_design', 'pmid']].drop_duplicates()

comb2 = pd.read_csv("results/full_superseries_with_pmid.csv", index_col=0)

comb2.columns = ['Unnamed: 0', 'series_id', 'characteristics_ch1',
       'extract_protocol_ch1', 'source_name_ch1', 'title_ex', 'title', 'summary',
       'overall_design', 'pmid']

super2 = comb2[['series_id', 'characteristics_ch1','title', 'summary',
       'overall_design', 'pmid']].drop_duplicates()

for_RAG = pd.concat([super1, super2])

for_RAG.to_csv("results/RAG.csv")