#This scripts converts aminoacid frequencies produced by workfolw/scripts/aa_frequency.py in a table of aminoacid frequencies for all species

import pandas as pd
import sys
#from collections import defaultdict

#Check if the expected number of arguments is provided
if len(sys.argv) < 4:
    print("Usage: python3 SCRIPT_NAME FILE_LIST INPUT_FOLDER OUTPUT_FOLDER")
    sys.exit(1)  # Exit the script with a non-zero status indicating error

#Open list of filenames (first argument of command line)
file_of_list_files = sys.argv[1]

#Get the second argument given in the command line and store it as input folder
input_folder = sys.argv[2]
#input_folder = '../../results/aa_frequencies/'

#Get the third argument given in the command line and store it as output folder
output_folder = sys.argv[3]

#Load list of files given in the command line
list_files = pd.read_csv(f'{file_of_list_files}', header = None, dtype=str)
list_files = list_files[0].tolist()

print("Processing files...")

#aa2freq = defaultdict(list) #cog2presence
aa2freq = {}

#Parsing irregularly-formatted cmscan file
for file_name in list_files:

    #Prepare to open file
    file_f = input_folder + str(file_name) + '.csv'
    df_tmp = pd.read_csv(file_f, sep=',', header=None)

    # Populate the dictionary
    for idx, row in df_tmp.iterrows():
        aminoacid = str(row[0])
        frequency = row[1]

        #Add info to dictionary -> matrix
        if(aminoacid in aa2freq):
            aa2freq[aminoacid][file_name] = frequency
        else:
            aa2freq[aminoacid] = {file_name: frequency}

#print(aa2freq)

#Convert to dataframe
aa_freq = pd.DataFrame.from_dict(aa2freq, orient='index')

#Substitute NA for zero
aa_freq = aa_freq.fillna(0)

#print(aa_freq)
print("Data saved in file", )

#Save the dataframe to a CSV file
aa_freq.to_csv(output_folder + 'aa_frequencies.csv', index=True)


