import pandas as pd

def check_if_superseries(df, gse_id):
    """check if the series is just the superseries
    """
    try:
        series_of_interest = df[df["Unnamed: 0"].isin([gse_id])]
        res = series_of_interest["summary"].str.contains("SuperSeries").bool()
        if res:
            return None
        else:
            return series_of_interest["title"].to_list()[0], series_of_interest["summary"].to_list()[0], series_of_interest["overall_design"].to_list()[0], series_of_interest["pubmed_id"].to_list()[0]
    except:
        return None

df_clean_stuff_remainder = pd.read_csv("results/combined_remainder.csv")
df_cleaned_superseries = df_clean_stuff_remainder.applymap(lambda x: remove_chars(x) if pd.notna(x) else x).drop_duplicates()

old_superseries = pd.read_csv("old_results/arch_geo_superseries.csv")
new_superseries = pd.read_csv("results/arch_geo_superseries.csv")
arch_remainder = new_superseries[~new_superseries["series_id"].isin(old_superseries["series_id"])]
df_for_annotation = arch_remainder

empty_entry1 = ("RNA-seq analysis after dBRD9-A treatment in human multiple myeloma cell line OPM2 and H929", \
               """BRD9 is a defining component of the non-canonical SWI/SNF complex, regulating gene expression by controlling 
               chromatin dynamics. We here identified high BRD9 expression as a poor prognostic factor in multiple myeloma (MM), 
               which is positively correlated with activation of ribosome biogenesis. Genetic and pharmacological depletion of BRD9 
               downregulates expression of ribosome biogenesis genes and disrupts protein synthesis maintenance machinery, thereby 
               inhibiting MM cell growth in both in vitro and in vivo preclinical models. Importantly, BRD9 interacts with BRD4 and 
               MYC to form a transcription initiation complex that promotes transcription of ribosomal biogenesis genes. 
               These results identify and validate BRD9 as a novel therapeutic target in MM which regulates ribosome biogenesis gene 
               transcription, and provide the framework for clinical evaluation of BRD9 degraders to improve patient outcome in MM. 
               This SuperSeries is composed of the SubSeries listed below.""", 
               "OPM2 and H929 cells were treated with either dBRD9-A or DMSO in triplicate for RNA-seq analysis.")


to_create = []

for i in tqdm(range(df_for_annotation.shape[0])):
    try:
        studies = df_for_annotation.iloc[i,1].split(",")
        original_list = [check_if_superseries(df_cleaned_superseries, i) for i in studies]
        filtered_list = [item for item in original_list if item is not None]
        to_create.append(filtered_list[0]) # tuple
    except:
        # print(traceback.print_exc())
        # print(studies)
        # print(df_for_annotation.iloc[i,1])
        if df_for_annotation.iloc[i,1] == 'GSE197487,GSE197492':
            to_create.append(empty_entry1)
        if df_for_annotation.iloc[i,1] == 'GSE197471,GSE197545':
            to_create.append(empty_entry2)
        if df_for_annotation.iloc[i,1] == 'GSE62829,GSE62830':
            to_create.append(empty_entry3)

superseries_added = pd.DataFrame(to_create)

superseries_added.columns = ["title", "summary", "overall_design", "pmid"]
res_nice = pd.concat([df_for_annotation.reset_index(), superseries_added], axis = 1, ignore_index=True)
res_nice.to_csv("results/superseries_remainder.csv")