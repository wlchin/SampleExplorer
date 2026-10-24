import argparse
from .embedding import download_model
from .sample_explorer import Query_DB
import os
import datetime
import sys
import logging
from importlib.metadata import PackageNotFoundError, version as package_version

def installed_version():
    """Installed package version, or 'unknown' when the metadata is missing (e.g. an uninstalled checkout)."""
    try:
        return package_version('sample-explorer')
    except PackageNotFoundError:
        return 'unknown'


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description='SampleExplorer Command Line Interface')
    parser.add_argument('--gene_list', type=str, nargs='?', default=None, help='Path to the gene list text file or comma-separated gene list (e.g. "IRF1,IRF2,IRF3")')
    parser.add_argument('--text_query', type=str, nargs='?', default=None, help='Path to a text file or a string, e.g. "studies with cardiac myocytes"')
    parser.add_argument('--archs4_file', type=str, help='Path to the ARCHS4 hdf5 file')
    parser.add_argument('--transcriptome_db', type=str, help='Path to the transcriptome database')
    parser.add_argument('--semantic_db', type=str, help='Path to the vector database')
    parser.add_argument('--search', type=str, help='Search strategy', default='semantic')
    parser.add_argument('--expand', type=str, help='Expansion strategy', default='transcriptome')
    parser.add_argument('--enrichment', action='store_true', default=False, help="Perform over-representation analysis (ORA) of each sample's top 5%% expressed genes")  # Added the --enrichment argument with default value False
    parser.add_argument('--usage', action='store_true', default=False, help='Print usage information')  # Added the --usage argument to print the usage information
    parser.add_argument('--download-model', action='store_true', default=False, help='Download the MiniLM language model once (e.g. before working offline) and exit')

    # Add any other command line arguments you need here

    parser.add_argument('--version', action='version', version='%(prog)s ' + installed_version())  # Added the --version argument to print the version

    if '--usage' in sys.argv:
        print("=====================================================")
        print("SampleExplorer Command Line Interface")
        print(f"Version {installed_version()}")
        print("Maintainer: WL Chin (melvin.chin@telethonkids.org.au)")
        print("=====================================================")
        print("")
        print("This is the usage information for the SampleExplorer Command Line Interface.")
        print("Please refer to the documentation on https://github.com/wlchin/SampleExplorer/ for detailed instructions on how to use this tool.")
        print("")
        print("Usage on command line:")
        print("----------------------")
        print("sample-explorer --gene_list IRF1,IRF2,IRF3 --text_query 'studies with cardiac myocytes' --archs4_file ARCHS4.h5 --transcriptome_db transcriptome.h5 --semantic_db vector.h5")
        print("")
        print("Usage using files:")
        print("------------------")
        print("sample-explorer --gene_list gene_list.txt --text_query text_query.txt --archs4_file ARCHS4.h5 --transcriptome_db transcriptome.h5 --semantic_db vector.h5")
        print("")
        print("Notes:")
        print("------")
        print("1. The gene names should be in uppercase")
        print("2. When using files, the gene_list.txt should contain a list of gene names (in UPPERCASE), one gene per line.")
        print("3. The file containing the text query should have the query string, which may contain multiple sentences, in plain text format")
        print("4. To perform enrichment (ORA), add the --enrichment flag")
        print("5. To download the language model once before first use (e.g. for offline machines), run: sample-explorer --download-model")
        print("")
        sys.exit(0)

    args = parser.parse_args()

    if args.download_model:
        for path in download_model():
            print(path)
        sys.exit(0)

    # Check if gene_list is provided
    if args.gene_list is not None:
        # Check if gene_list is a file or a comma-separated string
        if ',' in args.gene_list:
            # Split the comma-separated string into a list, dropping spaces and empty entries
            gene_list = [gene.strip() for gene in args.gene_list.split(',') if gene.strip()]
        else:
            # Read the gene list from the file line by line and clean each line
            # (utf-8-sig drops the byte order mark Windows editors may add)
            with open(args.gene_list, 'r', encoding='utf-8-sig') as gene_list_file:
                gene_list = [line.strip() for line in gene_list_file if line.strip()]
    else:
        gene_list = None

    # Check if text_query is provided
    if args.text_query is not None:
        # Check if text_query is a file or a string
        if args.text_query.endswith('.txt'):
            # Read the query from the file
            with open(args.text_query, 'r', encoding='utf-8-sig') as query_file:
                query = query_file.read()
        else:
            # Use the query as is
            query = args.text_query
    else:
        query = None

    new_query_db = Query_DB(args.semantic_db, args.transcriptome_db, args.archs4_file)
    result = new_query_db.search(geneset=gene_list, text_query=query, search=args.search, \
                                 expand=args.expand, perform_enrichment=args.enrichment)

    # Create the output folder if it doesn't exist
    output_folder = os.path.join('results', datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        logging.info(f"Results folder created: {output_folder}")

    # Save the result to CSV files in the output folder
    for table, filename in ((result.seed_studies, "relevant_seed_studies.csv"),
                            (result.expansion_studies, "relevant_expansion_studies.csv"),
                            (result.samples, "relevant_samples.csv")):
        if table is None:
            logging.info(f"No table to write for {filename}; skipping.")
            continue
        table.to_csv(os.path.join(output_folder, filename))
    logging.info(f"Results saved to {output_folder}")

if __name__ == '__main__':
    main()
