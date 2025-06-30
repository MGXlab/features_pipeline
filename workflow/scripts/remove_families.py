#This script filters cmscan reports. It keeps the original format, but filters out every line with hits that are non-prokaryotic Rfam families

#Prokaryotic Rfam families were found in: https://github.com/Rfam/rfam-taxonomy/blob/master/domains/bacteria.csv and https://github.com/Rfam/rfam-taxonomy/blob/master/domains/archaea.csv

#These two lists were downloaded on 30/06/2025 and stored in the Snakemake GITHub repository features_pipeline/config/archaea.csv and features_pipeline/config/bacteria.csv

#Usage:
#python3 script input_name output_name

import csv
import sys
import re

#Paths and names of Rfam family names
file_archaea = 'config/archaea.csv'
file_bacteria = 'config/bacteria.csv'

#Get list of Rfam families of archaea and bacteria from external files
with open(file_archaea, 'r') as f:
    reader = csv.reader(f)
    list_archaea = list({row[0].strip() for row in reader if row})

with open(file_bacteria, 'r') as f:
    reader = csv.reader(f)
    list_bacteria = list({row[0].strip() for row in reader if row})

#print(len(list_archaea))
#print(len(list_bacteria))

#Combine all Rfam family names
list_all = sorted(set(list_archaea + list_bacteria))
#print(list_all[:5])
#print(len(list_all))

#Get names of input and output files
if len(sys.argv) != 3:
    print("Usage: python3 script <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

#Open cmscan report and read it line by line
#Output the exact line, if the Rfam family matches an element of list_all created above
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        #Print comment lines to output file
        if line.startswith('#'):
            outfile.write(line)
        else:
            #Split lines of cmscan report by one or more spaces (the separator)
            fields = re.split(r'\s+', line.strip())
            #Print non-comment lines only if the Rfam family (3rd element of space separated line) matches any element of list_all
            if fields[2] in list_all:
                outfile.write(line)
            #else:
            #    print(line)
                                      
