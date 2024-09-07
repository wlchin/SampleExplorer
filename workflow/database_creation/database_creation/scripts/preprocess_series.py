import pandas as pd
from tqdm import tqdm

single_series = pd.read_csv("results/single_series.csv")
super_series = pd.read_csv("results/superseries.csv")

def remove_chars(input_string):
    """
    Removes specified characters from the input string.

    Args:
        input_string (str): The string from which characters will be removed.

    Returns:
        str: The input string with specified characters removed.
    """
    translation_table = str.maketrans('', '', "[]'")
    return input_string.translate(translation_table)

# Apply the function to every cell in the DataFrame
df_cleaned_singleseries = single_series.applymap(lambda x: remove_chars(x) if pd.notna(x) else x).drop_duplicates()
df_cleaned_superseries = super_series.applymap(lambda x: remove_chars(x) if pd.notna(x) else x).drop_duplicates()

def check_if_superseries(df, gse_id):
    """Check if the series is just the superseries.

    Args:
        df (pandas.DataFrame): The DataFrame containing the series information.
        gse_id (str): The GSE ID of the series to check.

    Returns:
        tuple or None: If the series is not a superseries, returns a tuple containing the title, summary, overall design, and PubMed ID of the series. 
                       If the series is a superseries, returns None.
    """
    series_of_interest = df[df["Unnamed: 0"].isin([gse_id])]
    res = series_of_interest["summary"].str.contains("SuperSeries").bool()
    if res:
        return None
    else:
        return series_of_interest["title"].to_list()[0], series_of_interest["summary"].to_list()[0], series_of_interest["overall_design"].to_list()[0], series_of_interest["pubmed_id"].to_list()[0]

def extract_singleseries(df, gse_id):
    """
    Extracts key points from the DataFrame based on the given GSE ID.

    Args:
        df (DataFrame): The DataFrame containing the data.
        gse_id (str): The GSE ID to filter the DataFrame.

    Returns:
        tuple: A tuple containing the title, summary, overall design, and PubMed ID of the series of interest.
    """
    series_of_interest = df[df["Unnamed: 0"].isin([gse_id])]
    res = series_of_interest["summary"].str.contains("SuperSeries").bool()
    return series_of_interest["title"].to_list()[0], series_of_interest["summary"].to_list()[0], series_of_interest["overall_design"].to_list()[0], series_of_interest["pubmed_id"].to_list()[0]

df_for_annotation = pd.read_csv("results/arch_geo_superseries.csv")
df_len = df_for_annotation.shape[0]
## these ones do not have series and cannot be retreived programmatically using GEOparse

#'GSE166953,GSE166955'
#'GSE197471,GSE197545'
#'GSE62829,GSE62830'

empty_entry1 = ("Withdrawn", "Withdrawn", "Withdrawn")
empty_entry2 = ('Cell division drives DNA methylation loss in late-replicating domains in primary human cells [RNA-seq]', 
               'We present the first experimental evidence that loss of DNA methylation at late-replicating regions, widely documented in aging and cancer, is directly driven by cell division. DNA hypomethylation, particularly at low-density CpGs in A:T-rich, partially methylated domains (PMD solo-WCGWs), tracks cumulative population doublings in primary cell culture. Cell cycle deceleration or full arrest resulted in a proportional decrease in the rate of DNA methylation loss. Lower culture oxygen conditions resulted in a slowed rate of methylation loss, suggesting that these CpGs may also be vulnerable to methylation loss during unscheduled, repair-associated methylation maintenance as well as scheduled, replication-associated methylation maintenance. By studying this methylation context in immortalized primary cells, we identified genomic and epigenomic features that are resistant or more susceptible to methylation loss. Finally, we leveraged this extensive DNA methylation dataset to develop a refined metric for the relative cumulative mitotic histories of human cells. These RNA-sequencing data accompany the larger DNA methylation dataset.',
               'RNA-sequencing was performed on total RNA extracted from serially-passaged human primary cells.')
empty_entry3 = ('Next-Generation Sequencing Analysis Reveals Differential Expression Profiles of miRNA-mRNA Target Pairs in KSHV-Infected Cells',
                """Purpose: Kaposi’s sarcoma (KS)-associated herpesvirus (KSHV) causes several lymphoproliferative disorders, including KS, a common AIDS-associated malignancy. Cellular and viral microRNAs (miRNAs) have been shown to play important roles in regulating the expression of genes in oncogenesis. Herpesviruses, including KSHV, encode for miRNAs that are involved in angiogenesis, inflammation and apoptosis. A better knowledge of the miRNA-mediated pathways that regulate KSHV infection is therefore essential for an improved understanding of viral infection and pathogenesis.
Methods: In this study, we used deep sequencing to analyze miRNA, both viral and human, and mRNA expression in KS tumor-derived human cells.
Results: This approach revealed 153 differentially expressed human miRNAs between KSHV-positive and -negative cells. Differential expression of eight miRNAs was independently confirmed by qRT-PCR. We additionally showed that a majority (~73%) of KSHV-regulated miRNAs are down-regulated, including most members of the 14q32 miRNA cluster. Specifically, human miR-409-3p, which is known to target the pro-angiogenic growth factor angiogenin and the inflammation marker fibrinogen-beta, was significantly down-regulated in KSHV-infected cells based on deep sequencing and qRT-PCR. Despite this substantial down-regulation of cellular miRNAs, hsa-miR-708-5p was significantly up-regulated by KSHV and has been shown to directly inhibit pro-apoptotic protease Caspase-2. Finally, we evaluated to what extent there was an inverse correlation between miRNA and mRNA expression levels. Using filtered datasets, we identified relevant canonical pathways that were significantly enriched.
Conclusion: Taken together, our data demonstrate that most human miRNAs affected by KSHV are repressed and our findings highlight the relevance of studying the post-transcriptional gene regulation of miRNAs for KSHV-associated malignancies.""",
               "Refer to individual Series. 6 samples analyzed (one cell type). Two experimental conditions: uninfected vs. chronically KSHV-infected cells (n=3). Two sequencing platforms: microRNA-Seq and mRNA-Seq.")


# merge annotations with the series

to_create = []

for i in tqdm(range(df_len)):
    try:
        studies = df_for_annotation.iloc[i,1].split(",")
        original_list = [check_if_superseries(df_cleaned_superseries, i) for i in studies]
        filtered_list = [item for item in original_list if item is not None]
        to_create.append(filtered_list[0]) # tuple
    except:
        if df_for_annotation.iloc[i,1] == 'GSE166953,GSE166955':
            to_create.append(empty_entry1)
        if df_for_annotation.iloc[i,1] == 'GSE197471,GSE197545':
            to_create.append(empty_entry2)
        if df_for_annotation.iloc[i,1] == 'GSE62829,GSE62830':
            to_create.append(empty_entry3)

superseries_added = pd.DataFrame(to_create)
superseries_added.columns = ["title", "summary", "overall_design", "pmid"]
pd.concat([df_for_annotation, superseries_added], axis = 1).to_csv("results/full_superseries_with_pmid.csv")

df_for_annotation = pd.read_csv("results/arch_geo_nosuperseries.csv")
df_len = df_for_annotation.shape[0]

to_create = []

for i in tqdm(range(df_len)):
    try:
        studies = df_for_annotation.iloc[i,1].split(",")
        original_list = [extract_singleseries(df_cleaned_singleseries, i) for i in studies]
        filtered_list = [item for item in original_list if item is not None]
        to_create.append(filtered_list[0]) # tuple
    except:
        studies = df_for_annotation.iloc[i,1]
        print(studies)

nosuperseries_added = pd.DataFrame(to_create)
nosuperseries_added.columns = ["title", "summary", "overall_design", "pmid"]

pd.concat([df_for_annotation, nosuperseries_added], axis = 1).to_csv("results/full_nosuperseries_with_pmid.csv")