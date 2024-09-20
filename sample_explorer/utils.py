import pandas as pd
import numpy as np
from .rag_embedding import Rag_embedding
import requests

class Sample_to_series_map:
    def __init__(self, h5_path, rag_index_path, rag_embedding_path):
        self.geo_metadata = self.process_metadata(h5_path)
        self.rag_embedding_object = Rag_embedding(rag_index_path, rag_embedding_path) 
        
    def process_metadata(self, h5_path):
        x = pd.read_pickle(h5_path)
        x["samples"] = x.index
        x = x[['series_id', 'samples']].drop_duplicates()
        return x
    
    def return_series_for_samples(self, samps):
        df = self.geo_metadata[self.geo_metadata["samples"].isin(samps)]
        return df
    
    def retrieve_text(self, gse_id_list):
        text_list = []
        for i in gse_id_list:
            try:
                txt = self.rag_embedding_object.get_text_linked_to_gse(i)
                text_list.append(txt)
            except:
                pass
        return text_list
    
    def calculate_similarites_from_samples(self, query, samples):
        samples = self.return_series_for_samples(samples)
        experimental_text = self.retrieve_text(samples["series_id"])
        avg = self.rag_embedding_object.get_averages_between_queries(query, experimental_text)
        return avg

    def calculate_similarites_from_series(self, query, series_list):
        experimental_text = self.retrieve_text(series_list)
        avg = self.rag_embedding_object.get_averages_between_queries(query, experimental_text)
        return avg


class ARCHS4_API_query:
    def __init__(self, h5_path):
        self.geo_metadata = self.process_metadata(h5_path)
        
    def process_metadata(self, h5_path):
        x = pd.read_pickle(h5_path)
        x["samples"] = x.index
        x = x[['series_id', 'samples']].drop_duplicates()
        return x
        
    def create_data_dict(self, direction, geneset):
        if direction == "up":
            data = {
                "type" : "geneset",
                "direction" : "similar",
                "signatureName" : "example_query",
                "species" : "human",
                "upgenes" : geneset
            }
            return data
        else:
            data = {
                "type" : "geneset",
                "direction" : "similar",
                "signatureName" : "example_query",
                "species" : "human",
                "downgenes" : geneset
            }
            return data
            
    def extract_sig(self, direction, sig_list):
        URL = 'https://maayanlab.cloud/rookpy/signature'
        data = self.create_data_dict(direction, sig_list)
        r = requests.post(URL, json=data)
        samps = r.json()["samples"]
        to_search = ["GSM" + str(i) for i in samps]
        return to_search
    
    def return_relevant_series_and_samples(self, direction, sig_list):
        samps = self.extract_sig(direction, sig_list)
        df = self.geo_metadata[self.geo_metadata["samples"].isin(samps)]
        return df, samps


class MsigDB_store:
    def __init__(self, signature_df_path, signature_metadata_df_path):
        self.signature_metadata_df = pd.read_csv(signature_metadata_df_path)
        self.signature_df = self.load_gene_set_database(signature_df_path)
    
    def get_gene_set(self, genesetname):
        index_loc = np.where(self.signature_df["gene_set_name"].isin([genesetname]))[0][0]
        genelist = self.signature_df.iloc[index_loc,2].split("\t")
        return genelist

    def get_title_from_index(self, genesetname):
        index_loc = np.where(self.signature_metadata_df["standard_name"].isin([genesetname]))[0][0]
        title_long = self.signature_metadata_df.iloc[index_loc,5]
        title_short = self.signature_metadata_df.iloc[index_loc,6]
        return title_short, title_long
    
    def load_gene_set_database(self, path):
        x = pd.read_csv(path, header = None)
        x.columns = ["unparsed_geneset"]
        correctly_formatted = x["unparsed_geneset"].str.split("\t", n = 2, expand = True)
        correctly_formatted.columns = ["gene_set_name", "link", "geneset"]
        return correctly_formatted
    
    def list_genesets(self):
        gene_set_names = self.signature_metadata_df["standard_name"]
        return gene_set_names
    
    def get_gene_set_by_name(self, geneset_name):
        gene_list = self.get_gene_set(geneset_name)
        short_title, long_title = self.get_title_from_index(geneset_name)
        response = {"geneset" : gene_list, "short_title" : short_title, "long_title" : long_title}
        return response
