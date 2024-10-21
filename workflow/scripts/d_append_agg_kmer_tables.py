#!/usr/bin/env python3

#Developer: Rose Brouns, modified by MB Walter Costa
'''Makes one table of all the kmer tables created with d_make_kmer_table.py
USAGE: SCRIPT input_file output_path'''
#input_file contains a list of IDs one per line and output path is a temporary folder

### HOUSEKEEPING

import pandas as pd
import sys
import os

### FUNCTION
def merge_kmer_profile(file_ids):
    '''
    Takes a list with the file ids as input and merge tables into dataframe.
    If a kmer is found duplicated, the merging went wrong 
    and the option to stop adding more tables is set.
    '''
    dfs = []  # List to hold dataframes before merging

    # Load and pre-process each table
    for file_name in file_ids:
        try:
            table = pd.read_csv(f'tmp/{file_name}_kmer{k}.tsv', sep='\t', index_col=0)
            table = table[~table.index.duplicated(keep='first')]
            dfs.append(table)
        except:
            print(f'>> {file_name} had an error')

    # Concatenate dataframes
    merged_df = pd.concat(dfs, axis=1)

    # Group by index and aggregate if needed
    merged_df = merged_df.groupby(merged_df.index).sum()  # Example: Sum values with same index

    return merged_df


if __name__ == '__main__':


    ### INPUT & OUTPUT

    # input_file is something like: list_kmer9_files.txt    
    input_file = sys.argv[1]
    k = input_file.split('_')[1].strip('kmer')
    #output_path should be output/
    output_path = sys.argv[2]
    
    ### RUN

    # Read the file_ids from the file
    with open(input_file, 'r') as file:
        file_ids = file.readlines()

    # Strip newline characters and remove empty lines
    file_ids = [line.strip() for line in file_ids if line.strip()]
    
    # Create empty DataFrame
    df = pd.DataFrame()

    df_merged = merge_kmer_profile(file_ids)
    
    # Check for duplicates
    if sum(df_merged.index.duplicated(keep=False)) > 0:
        print(f'> duplicates found in {chunk_name}')

    # Replace NaN by 0
    df_merged.fillna(0, inplace=True)
    df_merged = df_merged.apply(pd.to_numeric).astype(int)

    # Save to file with run-id
    df_merged.to_csv(f'{output_path}/kmer{k}_profiles.tsv', sep='\t')
