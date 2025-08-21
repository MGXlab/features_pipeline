#Script calculates aminoacid frequencies for the proteins of a protein FASTA file

from Bio import SeqIO
from collections import Counter
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Calculate amino acid frequencies from a FASTA file.")
parser.add_argument("fasta_file", help="Path to the FASTA file")
args = parser.parse_args()

# Initialize a counter to store total amino acid counts
aa_counter = Counter()

# Read all sequences from the FASTA file
with open(args.fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):

        #Remove stop signs '*' from the protein sequences, so that they do not get counted
        #print(record.seq)
        clean_sequence = str(record.seq).replace("*", "")
        #print(clean_sequence)

        aa_counter.update(str(clean_sequence))
        #break

# Total number of amino acids (excluding gaps or unknowns if needed)
total_aas = sum(aa_counter.values())

# Calculate frequencies
aa_frequencies = {aa: count / total_aas for aa, count in aa_counter.items()}

# Optional: print sorted results
for aa in sorted(aa_frequencies):
    print(f"{aa},{aa_frequencies[aa]:.4f}")

#print(f"Total frequency sum: {sum(aa_frequencies.values()):.4f}")
                                                       
