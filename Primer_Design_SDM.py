from pathlib import Path
from Bio import SeqUtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from primer3 import calc_tm # as calcTm - calcTm is deprecated

import pandas as pd
from itertools import combinations
import time
import os
import sys
import re
import math
from pathlib import Path
import argparse


# Design primers in 96-well format for Site-directed mutagenesis. One mutation at a time for this script.
# The primers are designed as 12 + 3 + 12. 
# If the location of mutagenesis site is too close to start or end, the length of primers can less than
# 27bp and in those cases manually verify the output. 
# change the k value in find_repeated_kmers(seq, k=16) to check repeats of different length
# Constants: Change as needed


### INPUTS ###########################
######################################

#CODON_TABLE_FILE = 'Primer_Design/Codon_Table_Standard.csv'
#ORF_FILE = 'Phytase_II.txt' # This file has intentional mistakes to test the script
#ORF_FILE = 'Primer_Design/HMT.txt'   ## Change input ORF filename
#MUTATION_LIST_FILE = 'Primer_Design/HMT_Plate2.csv' ## Change input mutation csv filename. column name should be "Mutations"


### OUTPUTS ###########################
######################################

#BASE_DIR = Path.cwd()
#print(BASE_DIR)
#OUTPUT_DIR = 'Primer_Design/Primers_HMT_Plate2' # Change output folder
#PRIMER_OUTPUT_FILE = 'HMT_Designed_primers.csv' # Change output file 1
#Forward_Primers_FILE = 'HMT_Forward_Primers_Plate2.csv' # Change output file 2
#Reverse_Primers_FILE = 'HMT_Reverse_Primers_Plate2.csv' # Change output file 3

# Create the full path to the file ###########################
######################################
def process_outputs(output_dir,primer_output_file,fwd_primers_file,rev_primers_file):
    base_dir=Path.cwd()
    Output_Path = base_dir / output_dir / primer_output_file
    Output_Path.parent.mkdir(parents=True, exist_ok=True)
    Output_Path_Fwd = base_dir / output_dir / fwd_primers_file
    Output_Path_Rev = base_dir / output_dir / rev_primers_file
    return Output_Path, Output_Path_Fwd, Output_Path_Rev

######################################
######################################

def remind_user_to_check_constants(mut_file,out_path,orf_file,codon_file):
    """Reminds the user to review constants and make changes if needed."""
    print("\n⚠️ Reminder: Please check the following inputs for correctness: \n")
    print(f"  - MUTATION_FILE: {mut_file}\n")
    print(f"  - PRIMER_OUTPUT: {out_path}\n")
    print(f"  - ORF_FILE: {orf_file}\n")
    print(f"  - CODON_TABLE: {codon_file}\n")
    print("Modify these and other relevant values in CLI arguments if necessary. Run with -h for help.\n")
        
        
def find_repeated_kmers(seq, k=16):
    """Check for repeated k-mers in the sequence. Change the k-value to check for 
    smaller repeats if needed."""
    seq_str = str(seq)  # Ensure we're working with a string
    kmers = [seq_str[i:i+k] for i in range(len(seq_str) - k + 1)]
    repeats = sum(seq_str.count(kmer) > 1 for kmer in kmers)
    
    if repeats:
        print(f'\nWarning: {repeats} repeats of {k}bp or more in the ORF. This will affect PCR & SDM.')

def read_orf_and_mutation_list(orf_file,mutation_list_file,codon_table_file):
    """Read the ORF sequence, mutation list, and codon table."""
    orf_seq = Seq(Path(orf_file).read_text().strip())
    mutation_list = pd.read_csv(mutation_list_file)
    codon_table = pd.read_csv(codon_table_file)
    return orf_seq, mutation_list, codon_table


def validate_position(amino_acid_pos, protein_len):
    """Ensure amino acid position is valid and within the ORF."""
    if not (1 <= amino_acid_pos <= protein_len):
        sys.exit('Amino acid position exceeds ORF bounds.')


def validate_mutations(mutation_list, orf_seq):
    """Check if the provided mutations align with translation of ORF sequence."""
    protein = orf_seq.translate()
    for mutation in mutation_list:
        if pd.notna(mutation) and mutation.strip() != "":
            pos = int(re.search(r'\d+', mutation).group()) - 1
            if protein[pos] != mutation[0]:
                sys.exit(f'ERROR: Mismatch at {pos + 1}. Expected {mutation[0]}, found {protein[pos]}.')
        else:
            sys.exit(f'ERROR: at null value in csv file. Format the csv file.')


def design_primers(orf_seq, mutations, codon_table):
    """Main primer design function."""
    primers = []

    for mutation in mutations['Mutations']:
        pos = int(re.search(r'\d+', mutation).group()) - 1
        validate_position(pos + 1, len(orf_seq) // 3)

        original_aa, target_aa = mutation[0], mutation[-1]
        codons = codon_table[codon_table['SingleLetter'] == target_aa]['Codon'].tolist()
        new_codon = find_optimal_codon(orf_seq, pos, codons)

        mutated_seq = orf_seq[:pos * 3] + new_codon + orf_seq[(pos + 1) * 3:]
        primer = extract_primer(mutated_seq, pos * 3)
        tm = int(math.ceil(calc_tm(str(primer), dv_conc=2, tm_method='santalucia', salt_corrections_method='owczarzy')))

        primers.append((mutation, primer, tm, int(SeqUtils.gc_fraction(primer)*100.0), len(primer))) #SeqUtils.GC is deprecated after Biopython 1.82 
        
    return pd.DataFrame(primers, columns=['Name', 'Sequence', 'Tm', 'GC', 'Length'])


def find_optimal_codon(orf_seq, pos, codons):
    """Find the most similar codon to the original one."""
    original_codon = str(orf_seq[pos * 3:(pos + 1) * 3])
    return max(codons, key=lambda codon: sum(a == b for a, b in zip(codon, original_codon)))


def extract_primer(seq, start):
    """Extract a 27bp primer around the mutation site."""
    return seq[max(0, start - 12):min(len(seq), start + 15)]


def create_primer_order_file(primers):
    """Create forward and reverse primers with appropriate names."""
    primer_order = []

    for _, row in primers.iterrows():
        fwd_name = f"Phy_{row['Name']}_Fwd"
        rev_name = f"Phy_{row['Name']}_Rev"
        rev_seq = str(Seq(row['Sequence']).reverse_complement())

        primer_order.extend([
            [fwd_name, row['Sequence'], row['Tm'], row['GC'], row['Length']],
            [rev_name, rev_seq, row['Tm'], row['GC'], row['Length']]
        ])

    return pd.DataFrame(primer_order, columns=['Name', 'Sequence', 'Tm', 'GC', 'Length'])


def separate_primers_by_type(primers, out_path_fwd, out_path_rev):
    """Separate forward and reverse primers into different CSV files."""

    # Create forward and reverse primer DataFrames safely
    fwd_primers = primers[primers['Name'].str.contains('_Fwd')].copy()
    rev_primers = primers[primers['Name'].str.contains('_Rev')].copy()

    # Define 96-well plate positions
    wells = [f"{row}{col}" for row in 'ABCDEFGH' for col in range(1, 13)]

    # Use .loc to assign wells to the DataFrames
    fwd_primers.loc[:, 'Well'] = wells[:len(fwd_primers)]
    rev_primers.loc[:, 'Well'] = wells[:len(rev_primers)]

    # Save to CSV files
    fwd_primers.to_csv(out_path_fwd, index=False)
    rev_primers.to_csv(out_path_rev, index=False)

    print("\nPrimers separated for IDT order in 96-well plate.")
    

def Validate_primer_length(primer_length):
    
    count = sum(1 for x in primer_length if x < 27)
    # If the count is greater than 0, print "Yes"
    if count:
        print(f'WARNING: {count} primers less than 27bp in the file.')

if __name__ == '__main__':
    start_time = time.time()

    parser = argparse.ArgumentParser()
    # MUTATION_LIST_FILE = 'Primer_Design/HMT_Plate2.csv'
    parser.add_argument('-m', '--Mutation_List', default='Primer_Design/HMT_Plate2.csv', help="List of Mutation Names mapped to well positions")
    # 'Primer_Design/Primers_HMT_Plate2'
    parser.add_argument('-o', '--Output_Directory', default='Primer_Design/Primers_HMT_Plate2', help="Output directory")
    # 'HMT_Designed_primers.csv'
    parser.add_argument('-f', '--Primer_Output_File', default='HMT_Designed_primers.csv', help="Output file for primers assoc. characteristic data")
    # 'HMT_Forward_Primers_Plate2.csv'
    parser.add_argument('-fwd', '--Forward_Primers_File', default='HMT_Forward_Primers_Plate2.csv', help="Output file for Forward primers")
    # Reverse_Primers_FILE = 'HMT_Reverse_Primers_Plate2.csv
    parser.add_argument('-rev', '--Reverse_Primers_File', default='HMT_Reverse_Primers_Plate2.csv', help="Output file for Reverse primers")

    # CODON_TABLE_FILE = 'Primer_Design/Codon_Table_Standard.csv'
    parser.add_argument('-c', '--Codon_Table_File', default='Primer_Design/Codon_Table_Standard.csv', help="Codon Translation Table Mapping File") 
    # ORF_FILE = 'Primer_Design/HMT.txt' 
    parser.add_argument('-orf', '--ORF_File', default='Primer_Design/HMT.txt', help="File containing Open Reading Frame Sequence") 

    args = parser.parse_args()

    out_path, path_fwd, path_rev = process_outputs(args.Output_Directory,args.Primer_Output_File,args.Forward_Primers_File,args.Reverse_Primers_File)

    remind_user_to_check_constants(args.Mutation_List,args.Output_Directory,args.ORF_File,args.Codon_Table_File)
    print(f'Working Directory: {os.getcwd()} \nProcessing...')
    
    orf_seq, mutations, codon_table = read_orf_and_mutation_list(args.ORF_File,args.Mutation_List,args.Codon_Table_File)

    find_repeated_kmers(orf_seq)
    #check if provided mutations align with translation
    validate_mutations(mutations['Mutations'].tolist(), orf_seq) 

    primers = design_primers(orf_seq, mutations, codon_table)
    Validate_primer_length(primers["Length"].tolist())
    
    primer_order = create_primer_order_file(primers)
    primer_order.to_csv(out_path, index=False)

    separate_primers_by_type(primer_order,path_fwd,path_rev)
    print(f"\nFinished in {time.time() - start_time:.6f} seconds.")
