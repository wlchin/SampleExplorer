# SampleExplorer [![codecov](https://codecov.io/gh/wlchin/bioRAG/graph/badge.svg?token=9GG94JA003)](https://codecov.io/gh/wlchin/bioRAG) ![workflow](https://github.com/wlchin/bioRAG/actions/workflows/python-package.yml/badge.svg) ![version](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11-blue) ![license](https://img.shields.io/badge/license-MIT-blue)

This is the repository for SampleExplorer. SampleExplorer can identify relevant studies within the ARCHS4 database, using a text-based query or a gene set. 

![Image Description](https://github.com/wlchin/bioRAG/blob/master/assets/BioRAG.png)

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [References](#references)

## Installation

SampleExplorer has been tested on python 3.9, 3.10, and 3.11. 

Install BioRAG via PyPI with the following command:

```
pip install sample_explorer
```

### Additional files


These additional files will need to be downloaded. Below is a table containing a description of these files, and the links.
| File Name | Description | Size | checksum |
|-----------|-------------|------|----------|
| semantic_db.h5ad   | [semantic vector store](https://data.pawsey.org.au/download/RNAseq_AB1_Renca/BioRAG/semantic_db.h5ad) | 58Mb | bd23f8835032dcd709f8ca05915038b3 |
| transcriptomic_db.h5ad   | [transcriptomic vector store](https://data.pawsey.org.au/download/RNAseq_AB1_Renca/BioRAG/transcriptomic_db.h5ad) | 38GB | d26e564653424ed6e189aea7d0be5d4d |
| human_gene_v2.2.h5    | ARCHS4[^1] [hdf5 database](https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.2.h5) | 37GB | f546fdecc0ba20bf2a5571628f132ca5 |

Both the semantic vector stores and the transcriptomic vector stores are AnnData[^2] objects, containing both the embedding matrices and the indices which link each embedding vector to an experiment in the ARCHS4 database. Additionally, the transcriptomic vector store also contains derived count data of "representative transcriptomes" in the ARCHS4 database. 

## Usage

To load SampleExplorer, follow these steps:

1. Download the additional files described above.
2. Instantiate a query_db object to perform SampleExplorer search in following way:
    
    ```python
    from sample_explorer.sample_explorer import Query_DB

    new_query_db = Query_DB(SEMANTIC_VECTOR_STORE_PATH, TRANSCRIPTOME_VECTOR_STORE_PATH, ARCHS4_HDF5_DATABASE_PATH)

    ```

3. Run SampleExplorer with a gene set query and/or textual query. The gene set is a python list, and the textual query is a python string.

    ```python

    text_query = "This text query usually describes the experiment"

    geneset_query = ["IFNG", "IRF1", "IFR2"]

    result = new_query_db.search(geneset = geneset_query, text_query = text_query)

    ```

4. The output is a Results object, which is composed of three pandas dataframes, which can be accessed via dot notation. 
- The "seed_studies" variable holds a dataframe containing study metadata from the search step. 
- The "expansion_studies" holds a dataframe of studies from the expansion step, 
- The "samples" variable holds a dataframe of the relevant samples, with metadata derived from the ARCHS4 database. 

```python

    result = new_query_db.search(geneset = geneset_query, text_query = text_query)

    # save the results as csv files
    result.seed_studies.to_csv("relevant_seed_studies.csv")
    result.expansion_studies.to_csv("relevant_expansion_studies.csv")
    result.samples.to_csv("relevant_samples.csv")

```

## Modifying searches using different inputs

All SampleExplorer searches should contain at least a text query or a gene set query. Users can choose one of several search strategies, with examples illustrated below:

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

By default, SampleExplorer will perform semantic search, followed by transcriptomic expansion. However, this search strategy can be modified. For instance, one can use semantic search, followed by perfoming the expansion step using semantic similiarity. In this case, the transcriptome vector store is not queried. 

The "seed" and "expand" parameters accepts either "transcriptome" or "semantic" as options.

If transcriptome search is performed as the initial step, the count data from "representative transcriptomes" in the transcriptome vector store is used to search for enriched samples, using ssGSEA (single sample gene set enrichment analysis) to rank the most relevant samples. In semantic search, SampleExplorer retrieves the most semantically similar studies to the text query, using cosine distance as the metric.

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

To further refine the set of samples and studies returned by SampleExplorer search, ssGSEA can be peformed on all samples returned by the query. Use the "perform_enrichment" parameter to specifiy if ssGSEA should be performed on all samples. If so, the returned dataframe will contain enrichment scores, pvalues and FDRs. The ssGSEA results will be stored as a dataframe in the "samples" attribute in the Results object.

```python

    result = new_query_db.search(text_query = text_query, geneset = None, search = "semantic", expand = "semantic", perform_enrichment = True)

    # save enrichment results from the sample dataframe
    results.samples.to_csv("samples_with_ssgsea_results.csv")

```


## License

SampleExplorer is published under the MIT License.

## Test datasets and continuous integration

The repository includes a [set](https://github.com/wlchin/SampleExplorer/tree/master/tests) of test data and testing scripts. The testing framework utilizes pytest.  

## Database creation and benchmarking workflows

The scripts for generating the embedding databases and performing benchmarking are located in the [workflow folder](./workflow/)

## References

[^1]: Lachmann A, Torre D, Keenan AB, Jagodnik KM, Lee HJ, Wang L, Silverstein MC, Ma'ayan A. Massive mining of publicly available RNA-seq data from human and mouse. Nat Commun. 2018 Apr 10;9(1):1366. doi: 10.1038/s41467-018-03751-6. PMID: 29636450; PMCID: PMC5893633.
[^2]: [AnnData](https://anndata.readthedocs.io/en/latest/)
