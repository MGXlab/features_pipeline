#!/bin/bash

# Check if the input file is provided in the script call
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <input> <folder>"
    exit 1
fi

# Input FASTA file
fasta_file="$1"

# Define name of output directory
output_dir="$2"

# Remove existing output directory if it exists
if [ -d "$output_dir" ]; then
    rm -r "$output_dir"
fi

# Create output directory 
mkdir -p "$output_dir"

# Define the number of proteins per file
proteins_per_file=500

# Count the total number of proteins in the input FASTA file
total_proteins=$(grep -c "^>" "$fasta_file")

# Calculate the number of files needed, considering the number of proteins per file
num_files=$(( (total_proteins + proteins_per_file - 1) / proteins_per_file ))

# Use awk to split the file and save each part to a separate file
awk -v proteins_per_file="$proteins_per_file" -v output_dir="$output_dir" -v num_files="$num_files" '
    BEGIN { protein_count = 0; file_count = 1; }
    /^>/ {
        if (protein_count % proteins_per_file == 0) {
            if (protein_count > 0) {
                close(output_file);
            }
            output_file = output_dir "/file" file_count ".faa";
            file_count++;
        }
        protein_count++;
    }
    { print >> output_file; }
    END { close(output_file); }
' "$fasta_file"

echo "$fasta_file file has been split into $num_files files in the $output_dir directory."

