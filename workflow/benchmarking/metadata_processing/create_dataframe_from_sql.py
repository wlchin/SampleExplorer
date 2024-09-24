import sqlite3
import pandas as pd

# Connect to the SQLite database file
connection = sqlite3.connect('data/msigdb_v2023.2.Hs.db')

# Create a cursor object to execute SQL queries
cursor = connection.cursor()

cursor.execute(
    """   SELECT collection_name, license_code, PMID AS PubMedID, GEO_id, description_full, description_brief, standard_name
    FROM gene_set gset
      INNER JOIN gene_set_details gsd ON gsd.gene_set_id = gset.id
      INNER JOIN publication pub ON pub.id = publication_id
    """
)

data = cursor.fetchall()
msigdb = pd.DataFrame(data)

columns = ["collection_name", "license_code", 
           "PMID", "GEO_id", "description_full", "description_brief", "standard_name"]

msigdb.columns = columns
msigdb.to_csv("data/gene_sets_msigdb.csv")