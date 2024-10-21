#!/usr/bin/env python3

#Developer: Rose Brouns, modified by MB Walter Costa
'''Makes from the fasta outputs with kmers counts from gerbil a table per sample
USAGE: SCRIPT input_file output_path k output_folder'''
#input_file contains a list of IDs one per line and output path is a temporary folder

### HOUSEKEEPING

import pandas as pd
import sys

### FUNCTION
def split_kmers_fasta(data):
    '''
    Put the kmers and kmer count from fasta format into 2 separate lists 
    '''
    kmers = []
    counts = []
    for line in data:
        if line.startswith('>'):
            # Take only number without > from line
            counts.append(int(line[1:]))
        else:
            kmers.append(line)
    
    return(kmers, counts)


if __name__ == '__main__':

    ### INPUT & OUTPUT

    # input_file is something like: list_kmer9_files.txt
    input_file = sys.argv[1]
    k =  sys.argv[3]
    output_folder = sys.argv[4]

    #output_path should be tmp/
    output_path = sys.argv[2]

    ### RUN

    # Read the file_ids from the file
    with open(input_file, 'r') as file:
        file_ids = file.readlines()

    # Strip newline characters and remove empty lines
    file_ids = [line.strip() for line in file_ids if line.strip()]
    
    #For each line/id extract kmers and counts
    for id in file_ids:

        # id is defined by Snakemake and has a strict format, e. g. CAADIQ000000000_kmer9.txt
        full_id = output_folder + '/kmer_files/' + id + '_kmer' + k + '.txt'
        
        # Read the data from the id file
        with open(full_id, 'r') as file:
            pre_data = file.read()

        # Split the data into lines and remove empty lines
        data = pre_data.split('\n')
        # Strip newline characters and remove empty lines        
        data = [i.strip() for i in data if i.strip()]

        # Put the kmers and counts in separate lists
        kmers, counts = split_kmers_fasta(data)
        
        # Create the DataFrame with kmers and counts
        df = pd.DataFrame({'': kmers, f'{id}': counts})

        #Save dataframe into temporary file
        df.to_csv(f'{output_path}/{id}_kmer{k}.tsv', sep='\t', index=False)
