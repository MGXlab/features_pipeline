#waltercostamb@gmail.com

#This script converts eggnog emapper's output to a binary CSV table. See wikis: 
# https://git.bia-christian.de/bia/lab_book_VEO/wiki/pipeline-of-features -> "Gene families"
# https://git.bia-christian.de/bia/lab_book_VEO/wiki/features -> "Gene functions"-> "Format eggonog output for SVM"

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

#Get the second argument given in the command line and store it as input folder
output_folder = sys.argv[3]

#Open list of filenames (first argument of command line)
file_of_list_files = sys.argv[1]

#Load list of files given in the command line
list_files = pd.read_csv(f'{file_of_list_files}', header = None, dtype=str)
list_files = list_files[0].tolist()

cog2presence = {}
line = 0

ts = time.time() 
print("Processing files...", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

#For each file, extract the 5th col. for every line/protein (eggNOG_OGs), create matrix presence/absence of orthologs
for file_name in list_files:
    
    #print(file_name)
    #Keep track of script, print processing of every 5th file
    #if(line %50 == 0):
    #    ts = time.time() 
    #    print("Processing file number", line, datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))
    
    #Open file, store 5th column
    t = input_folder + str(file_name) + '/' + str(file_name) + '.emapper.annotations'
    t2 = glob(t)
    #print(t)
    #print(t2)
    file1 = str(t2[0])

    #print(file1) 

    #Open file
    df = pd.read_csv(file1, sep = '\t', skiprows=5, header = None)
    #Store 5th column with the COGs (ortholog groups)
    #cog_list elements looks like: 
    #'COG4372@1|root,COG4372@2|Bacteria,4PM9B@976|Bacteroidetes,1J0UE@117747|Sphingobacteriia'
    
    cog_full_lst = df[4].tolist()
    
    #For each line/protein of the file, process it to store it in the presence/absence matrix
    for cog_full_el in cog_full_lst:
   
        #print(cog_full_el)

        #Avoid end of file lines begining with '#'
        if not isinstance(cog_full_el, float):
        
            #extract the highest level of COG (2nd element)
            cog_lst = cog_full_el.split('|')
        
            #remove extra string from 2nd element
            cog = cog_lst[1].replace('root,', '')
        
            id_file = file_name
            #print(cog_full_el,id_file)
        
            #Add info (ortholog's presence) to dictionary -> matrix
            if(cog in cog2presence):
                cog2presence[cog][id_file] = 1
            else:
                cog2presence[cog] = {id_file: 1}
                
    line = line + 1 
    
ts = time.time() 
print("Data in dictionary", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

#Convert to dataframe
df_cog = pd.DataFrame.from_dict(cog2presence, orient='index')

ts = time.time() 
print("Data in dataframe", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

#Substitute NA for zero
df_cog2 = df_cog.fillna(0)

#Convert data types from float to int
n = len(df_cog.columns)
df_cog = df_cog2.iloc[:, :n].astype(int)  

ts = time.time() 
print("Dataframe NA to 0", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

#Save the dataframe to a CSV file
df_cog.to_csv(output_folder + 'gene-family_profiles' + '.csv', index=True)

ts = time.time() 
print("Data saved in file", output_folder +"gene-family_profiles.csv", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))
