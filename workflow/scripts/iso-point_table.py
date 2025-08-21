#waltercostamb@gmail.com

#This script creates a CSV table of isoelectric points of ortholog proteins of different genomes.
#See wikis:
# https://git.bia-christian.de/bia/lab_book_VEO/wiki/pipeline-of-features -> "Isoelectric point of proteins"

#USAGE: (i) activate conda below to use python3 packages in draco and (ii) run the script
#conda activate /home/no58rok/tools/miniconda3/envs/bacterial_phenotypes
#python3 SCRIPT FILE_LIST COG_INPUT_FOLDER OUTPUT_FOLDER ISO_INPUT_FOLDER

#To create FILE_LIST, you can use the command lines below:
#$ls -lh INPUT_FOLDER/*annotations | sed 's/  */\t/g' | cut -f 9 | sed 's/\//\t/g' | cut -f3 > file_list.txt

import warnings
import re
import time
from datetime import datetime
from glob import glob
import pandas as pd
import numpy as np
import sys
import os

#Check if the expected number of arguments is provided
if len(sys.argv) < 5:
    print("Usage: python3 SCRIPT FILE_LIST COG_INPUT_FOLDER OUTPUT_FOLDER ISO_INPUT_FOLDER")
    sys.exit(1)  # Exit the script with a non-zero status indicating error

#Get argument from command line and store it as input folder
input_folder = sys.argv[2]

#Get argument from command line and store it as output folder
output_folder = sys.argv[3]

#Open list of filenames (first argument of command line)
file_of_list_files = sys.argv[1]

#Open argument from command line and store it as input folder contaning isoelectric points 
input_folder_iso = sys.argv[4]

#Load list of files given in the command line
list_files = pd.read_csv(f'{file_of_list_files}', header = None, dtype=str)
list_files = list_files[0].tolist()

print()
ts = time.time() 
print("Creating dictionary of protein names and isoelectric points...", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

#Create dictionary of protein names and isoelectric points
name2point = {}

#Go through list of files
for filename in list_files:

    #Form the complete path to the file
    file_path = os.path.join(input_folder_iso, f"{filename.split('/')[0]}-iso_point.csv")

    #Read the CSV file into a DataFrame
    iso_lines = pd.read_csv(file_path, header=0, sep=';', dtype=str)

    #print(iso_lines)

    for index,row in iso_lines.iterrows():
        prot_name = row['protein']
        iso_point = float(row['isoelectric_point'])

        #print(prot_name, iso_point)
        
        if filename not in name2point:
            name2point[filename] = {}

        #Populate dictionary
        if prot_name not in name2point[filename]:
            name2point[filename][prot_name] = iso_point

#print(name2point['1266999'])

ts = time.time() 
print("Processing files...", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

#Initialize dictionaries
cog2presence = {}
name2cog = {}

line = 0

#Create dictionary of cog2presence (isoelectric points of known genes, or, those that could be annotated by eggNOG emapper)
#For each file, extract the 5th col. for every line/protein (eggNOG_OGs), create matrix of orthologs, add isoelectric point of protein to fill the matrix
for file_name in list_files:
    
    #print(file_name)
    #Keep track of script, print processing of every 5th file
    #if(line %50 == 0):
    #    ts = time.time() 
    #    print("Processing file number", line, datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))
    
    #Open file, store 5th column (ortholog family name)
    t = input_folder + str(file_name) + '/' + str(file_name) + '.emapper.annotations'
    t2 = glob(t)
    #print(t)
    #print(t2)
    file1 = str(t2[0])

    #print(file1) 

    #Open file containing ortholog names
    df = pd.read_csv(file1, sep = '\t', skiprows=5, header = None)
    #Store 5th column with the COGs (ortholog groups)
    #cog_list elements looks like: 
    #'COG4372@1|root,COG4372@2|Bacteria,4PM9B@976|Bacteroidetes,1J0UE@117747|Sphingobacteriia'
    
    cog_full_lst = df[4].tolist()
    prot_name_lst = df[0]

    #For each line/protein of the file, process it to store isoelectric point matrix
    for i, cog_full_el in enumerate(cog_full_lst):
    #for cog_full_el in cog_full_lst:
   
        #print(cog_full_el)

        #Avoid end of file lines begining with '#'
        if not isinstance(cog_full_el, float):
        
            #extract prot_name
            prot_name = prot_name_lst[i]

            #extract the highest level of COG (2nd element)
            cog_lst = cog_full_el.split('|')
        
            #remove extra string from 2nd element
            pre_cog = cog_lst[1].replace('root,', '')
            #remove the NCBI tax ID
            cog = pre_cog.split('@')[0]

            #Populate dictionary with COG and gene name. This will be used in the following code block
            name2cog[prot_name] = cog

            id_file = file_name
            #print(cog_full_el,id_file)
        
            #Add isoelectric point to variable
            if (prot_name in name2point[id_file]):    
                iso_point = name2point[id_file][prot_name]
            else:
                iso_point = np.nan

            #Add info (isoelectric point) to dictionary -> matrix
            if(cog in cog2presence):
                cog2presence[cog][id_file] = iso_point
            else:
                cog2presence[cog] = {id_file: iso_point}
                
    line = line + 1 
   
ts = time.time() 
print("Data in dictionary", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

#Create a second dictionary, similar to cog2presence, but with all genes, not only those annotated by eggNOG emapper (or, those that contain a COG)
all2presence = cog2presence.copy()

#Further populate dictionary all2presence to contain isoelectric points of all genes, including those NOT annotated by eggNOG emapper
#First, go through every file_name
for file_name in list_files:

    #print(file_name)
    #print(name2point['1266999'])
    
    #Afterwards, go through every protein name of file $ID-iso_point.csv
    for prot_name, point in name2point[file_name].items():

        #Check if this name has a COG. If NOT, add element to dictionary
        if prot_name not in name2cog:

            iso_point = name2point[file_name][prot_name]
            #Add element
            all2presence[prot_name] = {file_name: iso_point} 

ts = time.time() 
print("Included proteins without COG annotation to second dictionary", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

#Convert to dataframe
df_cog = pd.DataFrame.from_dict(cog2presence, orient='index')
df_all = pd.DataFrame.from_dict(all2presence, orient='index')

ts = time.time() 
print("Data in dataframes", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

#Save the dataframes to two CSV files
df_cog.to_csv(output_folder + 'iso-points_profiles_known_orthologs' + '.csv', index=True)
df_all.to_csv(output_folder + 'iso-points_profiles_all_genes' + '.csv', index=True)

ts = time.time() 
print("Data 1 saved in file", output_folder +"iso-points_profiles_known_orthologs.csv", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))
print("Data 2 saved in file", output_folder +"iso-points_profiles_all_genes.csv", datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))
print()


