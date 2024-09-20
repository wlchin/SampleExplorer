#!/bin/bash

# Create a directory for downloads
download_dir="$(pwd)/downloads"
mkdir -p "$download_dir"
cd "$download_dir"

echo "Files will be downloaded to: $download_dir"

# Download files
files=(
    "https://data.pawsey.org.au/download/RNAseq_AB1_Renca/BioRAG/semantic_db.h5ad"
    "https://data.pawsey.org.au/download/RNAseq_AB1_Renca/BioRAG/transcriptomic_db.h5ad"
    "https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.2.h5"
)

for url in "${files[@]}"; do
    wget "$url"
    filename=$(basename "$url")
    echo "Downloaded $filename to: $download_dir/$filename"
done

# Install biorag
echo "Installing SampleExplorer..."
pip install sample-explorer

# Define expected md5sums
declare -A expected_md5sums=(
    ["semantic_db.h5ad"]="bd23f8835032dcd709f8ca05915038b3"
    ["transcriptomic_db.h5ad"]="d26e564653424ed6e189aea7d0be5d4d"
    ["human_gene_v2.2.h5"]="f546fdecc0ba20bf2a5571628f132ca5"
)

# Check md5sums
echo "Checking MD5 checksums..."
for file in "${!expected_md5sums[@]}"; do
    actual_md5=$(md5sum "$file" | awk '{print $1}')
    expected_md5=${expected_md5sums[$file]}
    
    if [ "$actual_md5" = "$expected_md5" ]; then
        echo "$file: MD5 checksum matches"
    else
        echo "$file: MD5 checksum mismatch!"
        echo "  Expected: $expected_md5"
        echo "  Actual:   $actual_md5"
    fi
done

echo "All operations completed. Files are stored in: $download_dir"