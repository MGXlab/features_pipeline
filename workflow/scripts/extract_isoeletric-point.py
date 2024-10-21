import re
import sys

#This script extracts the peptide/protein name and isoeletric point of an output file of EMBOSS pepstats

# Check if a file name is provided as a command-line argument, otherwise report warning
if len(sys.argv) != 2:
    print("Usage: python3 script.py <input>")
    sys.exit(1)

# Get the file name from the command-line argument
file_name = sys.argv[1]

# Read the content of the file
try:
    with open(file_name, 'r') as file:
        input_text = file.read()
except FileNotFoundError:
    print(f"Error: File '{file_name}' not found.")
    sys.exit(1)

# Define a regular expression pattern to match protein names and Isoelectric Point
pep = re.compile(r'PEPSTATS of (\S+) from')
ie = re.compile(r'Isoelectric Point = (\d+\.\d+)')

# Find all matches in the input text
matches_pep = pep.findall(input_text)
matches_ie = ie.findall(input_text)

print("protein;isoelectric_point")

# Print extracted information
for i in range(len(matches_pep)):
    print(f"{matches_pep[i]};{matches_ie[i]}")
