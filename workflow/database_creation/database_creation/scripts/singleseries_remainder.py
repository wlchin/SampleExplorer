import pandas as pd

def extract_singleseries(df, gse_id):
    """return key points
    """
    series_of_interest = df[df["Unnamed: 0"].isin([gse_id])]
    res = series_of_interest["summary"].str.contains("SuperSeries").bool()
    return series_of_interest["title"].to_list()[0], series_of_interest["summary"].to_list()[0], series_of_interest["overall_design"].to_list()[0], series_of_interest["pubmed_id"].to_list()[0]

df_clean_stuff_remainder = pd.read_csv("./combined_remainder.csv")
df_cleaned_superseries = df_clean_stuff_remainder.applymap(lambda x: remove_chars(x) if pd.notna(x) else x).drop_duplicates()

old_superseries = pd.read_csv("../old_results/arch_geo_nosuperseries.csv")
new_superseries = pd.read_csv("../results/arch_geo_nosuperseries.csv")
arch_remainder = new_superseries[~new_superseries["series_id"].isin(old_superseries["series_id"])]
df_for_annotation = arch_remainder

empty_entry1 = ("Withdrawn", "Withdrawn", "Withdrawn")

to_create = []

for i in tqdm(range(df_for_annotation.shape[0])):
    try:
        studies = df_for_annotation.iloc[i,1].split(",")
        original_list = [extract_singleseries(df_cleaned_superseries, i) for i in studies]
        filtered_list = [item for item in original_list if item is not None]
        to_create.append(filtered_list[0]) # tuple
    except:
        # print(traceback.print_exc())
        # print(df_for_annotation.iloc[i,1])
        if df_for_annotation.iloc[i,1] == 'GSE222563':
            to_create.append(empty_entry1)
        # if df_for_annotation.iloc[i,1] == 'GSE197471,GSE197545':
        #     to_create.append(empty_entry2)
        # if df_for_annotation.iloc[i,1] == 'GSE62829,GSE62830':
        #     to_create.append(empty_entry3)

superseries_added = pd.DataFrame(to_create)

superseries_added.columns = ["title", "summary", "overall_design", "pmid"]
res_nice = pd.concat([df_for_annotation.reset_index(), superseries_added], axis = 1, ignore_index=True)

res_nice.to_csv("singleseries_remainder.csv")