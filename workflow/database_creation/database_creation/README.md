
## Snakemake workflow for BioRAG databases

### Introduction

This repository contains the snakemake workflow for generating the embedding databases for BioRAG. The python package BioRAG is available at this [link](https://github.com/wlchin/bioRAG).

BioRAG requires two vector databases to conduct retrival of relevant studies. The snakemake workflow details the creation of these databases from human RNAseq data in the ARCHS4 hdf5 files.

### Software requirements.

The following packages are required:
1. archs4py - to retrieve data from hdf5 ARCHS4 files
2. geoparse - to retrieve GEO metadata from ARCHS4 studies
3. Sentence-transformers - to create vector embeddings from GEO metadata
4. scikit-learn - to create transcriptomic embeddings from hdf5 ARCHS4 files
4. numpy - for persistence of transformed, memory-mapped count data
5. pandas - for data preproceesing of GEO metadata
6. snakemakem - for worfklow construction.

### Running the workflow

Please refer to [documentation](https://snakemake.readthedocs.io/en/stable/) to run the workflow from scratch. The workflow was created using snakemake version 6.0