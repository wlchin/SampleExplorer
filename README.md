# BioRAG

This is the repository for BioRAG. BioRAG can identify gene signature enrichment under relevant transcriptomic conditions within the ARCHS4 database, using a researcher-supplied gene set or text-based experimental description. 

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [License](#license)

## Installation

Install BioRAG via PyPI with the following command:

    ```python
    pip install biorag
    ```

### Additional files

These additional files will need to be downloaded. Below is a table containing a description of these files, and the links.

| File Name | Description | Link |
|-----------|-------------|------|
| semantic_db    | Description 1 | [Link 1](https://example.com/file1) |
| transcriptomic_db   | Description 2 | [Link 2](https://example.com/file2) |
| ARCHS4 count data    | Description 3 | [Link 3](https://example.com/file3) |

## Usage

Examples and instructions on how to use the project.

To load BioRAG, follow these steps:

1. Download the following files.
2. Initialize the BioRAG object in the following way:
3. Run biorag with a gene set query and textual query:

    ```bash
    biorag
    ```

The gene set is a list, and the textual query is a string python object.

## License

BioRAG is provided under the GNU General Public License (GPL), which ensures that users have the freedom to run, study, modify, and distribute the software.