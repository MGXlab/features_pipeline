#waltercostamb@gmail.com

#This script converts Infernal's cmscan output to a binary CSV table. See wiki: 
# https://git.bia-christian.de/bia/lab_book_VEO/wiki/pipeline-of-features-rules -> "Riboswitches"

#USAGE: (i) activate conda below to use python3 packages in draco and (ii) run the script
#conda activate bacterial_phenotypes
#python3 SCRIPT FILE_LIST INPUT_FOLDER OUTPUT_FOLDER

#To create FILE_LIST, you can use the command lines below:
#$ls -lh INPUT_FOLDER/*annotations > pre_file_list.txt
#$sed 's/  */\t/g' pre_file_list.txt | cut -f 9 | sed 's/\//\t/g' | cut -f3 > file_list.txt

import warnings
import re
import time
from datetime import datetime
from glob import glob
import pandas as pd
import sys

#Check if the expected number of arguments is provided
if len(sys.argv) < 4:
    print("Usage: python3 genes_table.py FILE_LIST INPUT_FOLDER OUTPUT_FOLDER")
    sys.exit(1)  # Exit the script with a non-zero status indicating error

#Get the second argument given in the command line and store it as input folder
input_folder = sys.argv[2]
#input_folder = '../../results/ncRNAs_infernal/'

#Get the second argument given in the command line and store it as input folder
output_folder = sys.argv[3]

#Open list of filenames (first argument of command line)
file_of_list_files = sys.argv[1]

#Load list of files given in the command line
list_files = pd.read_csv(f'{file_of_list_files}', header = None, dtype=str)
list_files = list_files[0].tolist()
#with open('../../config/files.txt', 'r') as file:
#    list_files = file.read().splitlines()
#list_files

print("Processing files...")

cog2presence = {}
cog2presence_binary = {}

#Parsing irregularly-formatted cmscan file
for file_name in list_files:
    
    #Prepare to open file
    t = input_folder + str(file_name) + '.cmscan'
    with open(t, 'r') as file:
        # Read all file lines and store to list
        lines = file.readlines()
        # strip newline characters of lines
        lines = [line.strip() for line in lines]

    #Marker set to 0 unless it finds a line to skip
    marker = 0

    #Go through every line to properly parse the lines
    for line in lines:
        
        #If header is found, get the following column names, there should be 29 columns
        if 'idx' in line:

            colnames = ['idx' ,'target name', 'accession1', 'query name', 'accession2', 'clan name', 'mdl', 'mdl from', 'mdl to', 'seq from',
                       'seq to', 'strand', 'trunc', 'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'olp', 
                       'anyidx', 'afrct1', 'afrct2', 'winidx', 'wfrct1', 'wfrct2', 'mdl len', 'seq len', 'description of target']

            #Start df
            df = pd.DataFrame(columns=colnames)
        
        #If line starts with # skip it (except for the header, which has already been parsed above)
        elif line.startswith('#'):
            marker = 1
        else:
            # Split the line by one or more spaces (\s+), get first 28 elements
            col_elements = re.split(r'\s+', line.strip())[:28]  # strip() removes any trailing/leading whitespaces
            
            #For 29th element get all elements and join them by space
            pre_description = re.split(r'\s+', line.strip())[28:]
            description = ' '.join(pre_description)
            col_elements.append(description)
            
            # Add the first 29 elements as a new row in the DataFrame
            df.loc[len(df)] = col_elements
            # Drop un-wanted columns
            df_final = df.copy()
            df_final = df_final.drop(columns=['accession2', 'clan name', 'mdl', 'mdl from', 'mdl to', 
                                   'pass', 'score', 'inc', 'olp', 
                                   'anyidx', 'afrct1', 'afrct2', 'winidx', 'wfrct1', 'wfrct2', 'mdl len'])
    
    #Filter for false positive hits of cmscan
    df_final['bias'] = pd.to_numeric(df_final['bias'])
    df_final = df_final[df_final['bias'] <= 50]
    df_final['E-value'] = pd.to_numeric(df_final['E-value'])
    df_final = df_final[df_final['E-value'] < 0.001]   
    
    #df_final.head()
    #df_final.shape

    # Get the value counts for the 'target name'/ncRNA column and convert it to a dictionary
    target_name_counts = df_final['target name'].value_counts().to_dict()
    #target_name_counts
    
    #For each ncRNA of the file, process it to store it in the presence/absence or count matrix cog2presence
    for ncRNA in target_name_counts:
           
        id_file = file_name
        
        #Add info (ortholog's presence) to dictionary -> matrix
        if(ncRNA in cog2presence):
            cog2presence[ncRNA][id_file] = target_name_counts[ncRNA]
            cog2presence_binary[ncRNA][id_file] = 1
        else:
            cog2presence[ncRNA] = {id_file: target_name_counts[ncRNA]}
            cog2presence_binary[ncRNA] = {id_file: 1}

#cog2presence

print("Data in dictionary")

#Convert to dataframe
df_cog = pd.DataFrame.from_dict(cog2presence, orient='index')
df_cog_binary = pd.DataFrame.from_dict(cog2presence_binary, orient='index')

print("Data in dataframe")

#Substitute NA for zero
df_cog2 = df_cog.fillna(0)
df_cog2_binary = df_cog_binary.fillna(0)

#Convert data types from float to int
n = len(df_cog.columns)
df_cog = df_cog2.iloc[:, :n].astype(int)  

n_binary = len(df_cog_binary.columns)
df_cog_binary = df_cog2_binary.iloc[:, :n_binary].astype(int)  

print("Dataframe NA to 0")

#df_cog

#Save the dataframe to a CSV file
df_cog.to_csv(output_folder + 'ncRNA_profiles_counts' + '.csv', index=True)
df_cog_binary.to_csv(output_folder + 'ncRNA_profiles_binary' + '.csv', index=True)

print("Data saved in file", output_folder +"ncRNA_profiles_counts.csv")
print("Data saved in file", output_folder +"ncRNA_profiles_binary.csv")
