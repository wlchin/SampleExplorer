# BioRAG

This is the repository for BioRAG. BioRAG can identify relevant studies within the ARCHS4 database, using a text-based query or a gene set. 

![Image Description](https://github.com/wlchin/bioRAG/blob/master/assets/BioRAG.png)

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [License](#license)

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
| human_gene_v2.2.h5    | ARCHS4 count data | [Link 3](https://example.com/file3) |

## Usage

To load BioRAG, follow these steps:

1. Download the additional files described above.
2. Instantiate a query_db object to perform BioRAG search in following way:
    
    ```python
    from biorag import query_db

    new_query_db = query_db(SEMANTIC_VECTOR_STORE_PATH, TRANSCRIPTOME_VECTOR_STORE_PATH, COUNT_H5_PATH)

    ```

3. Run BioRAG with a gene set query and textual query. The gene set is a python list, and the textual query is a python string.

    ```python

    text_query = "This text query usually describes the experiment"

    geneset_query = ["IFNG", "IRF1", "IFR2"]

    result = new_query_db.search(geneset = geneset_query, text_query = text_query)

    ```

4. The results from BioRAG are in the form of a a pandas dataframe. The rows are <u>samples</u> from the most relevant studies to the query with associated experimental metadata.

## Modifying searches using different inputs

All BioRAG searches should contain at least a text query. Users can choose one of several search strategies:

- Strategy 1: Text as input only

    ```python

    result = new_query_db.search(text_query = text_query, geneset = None)

    ```

- Strategy 2: Text and gene set as input

    ```python

    result = new_query_db.search(text_query = text_query, geneset = geneset_query)

    ```

- Strategy 3: Gene set as input

    ```python

    result = new_query_db.search(text_query = None, geneset = geneset_query)

    ```

## Modifying searches using different expansion strategies

By default, BioRAG will perform semantic search, followed by transcriptomic expansion. However, this search strategy can be modified to use only semantic search alone (perfoming the expansion step using semantic similiarity). The "seed" and "expand" parameters accepts either "transcriptome" or "semantic" as options.

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

The types of searches performed depends on input. If only a gene set is supplied, the search step defaults to "transcriptome", with the user able to select between "transcriptome" and "semantic" for the expansion step.

## Optional single sample GSEA (ssGSEA)

To further refine the set of samples and studies returned by BioRAG search, ssGSEA can be peformed on all samples returned by the query. Use the "perform_enrichment" parameter to specifiy if ssGSEA should be performed on all samples. If so, the returned dataframe will contain enrichment scores, pvalues and FDRs.

```python

    result = new_query_db.search(text_query = text_query, geneset = None, search = "semantic", expand = "semantic", perform_enrichment = True)

```

## License

BioRAG is provided under the GNU General Public License (GPL).