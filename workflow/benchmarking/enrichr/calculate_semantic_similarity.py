from biorag.biorag import Rag_embedding
from biorag.utils import MsigDB_store
import pandas as pd
from tqdm import tqdm
import traceback
import os

msigdb_store = MsigDB_store("data/msigdb.v2023.2.Hs.symbols.gmt", "data/test_set_misgdb.csv")

rag = Rag_embedding("data/important_data_for_LLM/rag_index_v2.pkl", 
                    "data/important_data_for_LLM/rag_embedding_matrix_v2.pkl")

genesets = os.listdir("tempfiles/results_enrichr_metadata/") 

# remove the extension name from the strings in the list
filtered_list = [i.replace(".csv", "") for i in genesets]

dicto = []

for i in tqdm(filtered_list):
    try:
        print(i)
        testset = msigdb_store.get_gene_set_by_name(i)
        query = testset["long_title"]
        res_dict = {}
        res_dict["signature"] = i
        filename_source = "tempfiles/results_enrichr_metadata/" + i + ".csv"
        x = pd.read_csv(filename_source)
        text_list = []
        num_rows = x.shape[0]
        for i in range(num_rows):
            text = eval(x.iloc[i,1])[0] + eval(x.iloc[i,7])[0] + eval(x.iloc[i,8])[0] # this draws from the metadata columns
            text_list.append(text)
        av = rag.get_averages_between_queries(query, text_list)
        res_dict["average_sim"] = av
        dicto.append(res_dict)
    except:
        print("failure", "batch: ", res_dict["signature"])
        print(traceback.print_exc())
    
pd.DataFrame(dicto).to_csv("tempfiles/cosine_similarity_enrichr_queries.csv")