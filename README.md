# BioRAG

This is the repository for BioRAG[^1]. BioRAG can identify relevant studies within the ARCHS4 database, using a text-based query or a gene set. 

![Image Description](https://github.com/wlchin/bioRAG/blob/master/assets/BioRAG.png)

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [References](#references)

## Installation

Install BioRAG via PyPI with the following command:

```
pip install biorag
```

### Additional files

These additional files will need to be downloaded. Below is a table containing a description of these files, and the links.

| File Name | Description | Link |
|-----------|-------------|------|
| semantic_db.h5ad   | semantic vector store | [Link 1](https://example.com/file1) |
| transcriptomic_db.h5ad   | transcriptomic vector store | [Link 2](https://example.com/file2) |
| human_gene_v2.2.h5    | ARCHS4 hdf5 database[^2] | [Link 3](https://example.com/file3) |

Both the semantic vector stores and the transcriptomic vector stores are AnnData objects, containing both the embedding matrices and the indices which link each embedding vector to an experiment in the ARCHS4 database. Additionally, the transcriptomic vector store also contains derived count data of "representative transcriptomes" in the ARCHS4 database. If transcriptome search is performed as the initial step, the count data is used to search for studies in the ARCHS4 database containing enriched samples, using ssGSEA (single sample gene set enrichment analysis) to rank the most relevant samples. 

## Usage

To load BioRAG, follow these steps:

1. Download the additional files described above.
2. Instantiate a query_db object to perform BioRAG search in following way:
    
    ```python
    from biorag.biorag import Query_DB

    new_query_db = Query_DB(SEMANTIC_VECTOR_STORE_PATH, TRANSCRIPTOME_VECTOR_STORE_PATH, ARCHS4_HDF5_DATABASE_PATH)

    ```

3. Run BioRAG with a gene set query and/or textual query. The gene set is a python list, and the textual query is a python string.

    ```python

    text_query = "This text query usually describes the experiment"

    geneset_query = ["IFNG", "IRF1", "IFR2"]

    result = new_query_db.search(geneset = geneset_query, text_query = text_query)

    ```

4. The output is a Results object, which is composed of three pandas dataframes, which can be accessed via dot notation. The "seed_studies" variable holds a dataframe containing study metadata from the search step. The "expansion_studies" holds a dataframe of studies from the expansion step, whilst the "samples" variable holds a dataframe of the relevant samples, with metadata derived from the ARCHS4 database. Below, these dataframes are accessed from the Results object and saved as CSV files.

    ```python

    result = new_query_db.search(geneset = geneset_query, text_query = text_query)

    result.seed_studies.to_csv("relevant_seed_studies.csv")
    result.expansion_studies.to_csv("relevant_expansion_studies.csv")
    result.samples.to_csv("relevant_samples.csv")

    ```

## Modifying searches using different inputs

All BioRAG searches should contain at least a text query or a gene set query. Users can choose one of several search strategies, with examples illustrated below:

- Strategy 1: Text only as input

    ```python

    result = new_query_db.search(text_query = text_query, geneset = None)

    ```

- Strategy 2: Text and gene set as input

    ```python

    result = new_query_db.search(text_query = text_query, geneset = geneset_query)

    ```

- Strategy 3: Gene set only as input

    ```python

    result = new_query_db.search(text_query = None, geneset = geneset_query)

    ```

## Modifying searches using different expansion strategies

By default, BioRAG will perform semantic search, followed by transcriptomic expansion. However, this search strategy can be modified. For instance, one can use semantic search, followed by perfoming the expansion step using semantic similiarity. In this case, the transcriptome vector store is not queried. 

The "seed" and "expand" parameters accepts either "transcriptome" or "semantic" as options.

1. Example 1 - perform semantic search, followed by transcriptome expansion.

    ```python

    result = new_query_db.search(text_query = text_query, geneset = geneset_query, search = "semantic", expand = "transcriptome")

    ```

2. Example 2 - perform transcriptome search, followed by semantic expansion.

    ```python

    result = new_query_db.search(text_query = text_query, geneset = geneset_query, search = "transcriptome", expand = "semantic")

    ```

2. Example 3 - perform semantic search, followed by semantic expansion, using only a text query.

    ```python

    result = new_query_db.search(text_query = text_query, geneset = None, search = "semantic", expand = "semantic")

    ```

The types of searches which can be performed depend on input. For instance,if only a gene set is supplied, the search step defaults to "transcriptome", with the user able to select between "transcriptome" and "semantic" for the expansion step.

## Optional single sample gene set enrichment analysis (ssGSEA)

To further refine the set of samples and studies returned by BioRAG search, ssGSEA can be peformed on all samples returned by the query. Use the "perform_enrichment" parameter to specifiy if ssGSEA should be performed on all samples. If so, the returned dataframe will contain enrichment scores, pvalues and FDRs. The ssGSEA results will be stored as a dataframe in the "samples" attribute in the Results object.

```python

    result = new_query_db.search(text_query = text_query, geneset = None, search = "semantic", expand = "semantic", perform_enrichment = True)

    results.samples.to_csv("samples_with_ssgsea_results.csv")

```

## License

Creative Commons Attribution 4.0 International. The Creative Commons Attribution license allows re-distribution and re-use of a licensed work on the condition that the creator is appropriately credited. 

## References

[^1]: Chin WL, & Lassmann, T. (2024). Language models improve the discovery of public RNA-seq data. *bioRxiv*.
[^2]: Lachmann A, Torre D, Keenan AB, Jagodnik KM, Lee HJ, Wang L, Silverstein MC, Ma'ayan A. Massive mining of publicly available RNA-seq data from human and mouse. Nat Commun. 2018 Apr 10;9(1):1366. doi: 10.1038/s41467-018-03751-6. PMID: 29636450; PMCID: PMC5893633.
