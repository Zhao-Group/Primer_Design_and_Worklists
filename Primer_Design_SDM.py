from pathlib import Path
from Bio import SeqUtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from primer3 import calc_tm # as calcTm - calcTm is deprecated
import primer3

from Bio.SeqUtils import MeltingTemp as mt

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

def remind_user_to_check_constants(mut_file,out_path,orf_file,codon_file_id):
    """Reminds the user to review constants and make changes if needed."""
    print("\n⚠️ Reminder: Please check the following inputs for correctness: \n")
    print(f"  - MUTATION_FILE: {mut_file}\n")
    print(f"  - PRIMER_OUTPUT: {out_path}\n")
    print(f"  - ORF_FILE: {orf_file}\n")
    print(f"  - CODON_TABLE: {CodonTable.unambiguous_dna_by_id[int(codon_file_id)]}\n")
    print("Modify these and other relevant values in CLI arguments if necessary. Run with -h for help.\n")
        
        
def find_repeated_kmers(seq, k=16):
    """Check for repeated k-mers in the sequence. Change the k-value to check for 
    smaller repeats if needed."""
    seq_str = str(seq)  # Ensure we're working with a string
    kmers = [seq_str[i:i+k] for i in range(len(seq_str) - k + 1)]
    repeats = sum(seq_str.count(kmer) > 1 for kmer in kmers)
    
    if repeats:
        print(f'\nWarning: {repeats} repeats of {k}bp or more in the ORF. This will affect PCR & SDM.')

def read_orf_and_mutation_list(orf_file,mutation_list_file,codon_table_file_id):
    """Read the ORF sequence, mutation list, and codon table."""
    orf_seq = Seq(Path(orf_file).read_text().strip())

    # ✅ Check if cDNA starts with start codon ATG
    if not orf_seq.upper().startswith("ATG"):
        print("⚠️ Warning: The provided cDNA sequence does not start with ATG (start codon).")
        
    mutation_list = pd.read_csv(mutation_list_file)
    #codon_table = CodonTable.unambiguous_dna_by_id[codon_table_file_id].tolist() #pd.read_csv(codon_table_file)
    return orf_seq, mutation_list #, codon_table


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


def design_primers(orf_seq, mutations, codon_table_id):
    """Main primer design function."""
    primers = []

    for mutation in mutations['Mutations']:
        pos = int(re.search(r'\d+', mutation).group()) - 1
        validate_position(pos + 1, len(orf_seq) // 3)

        # Get all possible translation table combinations for the aa resiudes in question
        original_aa, target_aa = mutation[0], mutation[-1]
        codon_table = CodonTable.unambiguous_dna_by_id[int(codon_table_id)].forward_table
        codons = [key for key, value in codon_table.items() if value == target_aa]
        #codons = codon_table[codon_table['SingleLetter'] == target_aa]['Codon'].tolist()
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


def get_Amino_Acid_Sequence(seq):
    """Translate a DNA sequence to its corresponding amino acid sequence."""

    'Biopython - https://biopython.org/docs/1.75/api/Bio.Seq.html'
    return str(seq.translate())

def get_Stop_Codon(seq):
    """Check if the DNA sequence contains a stop codon."""
    seq=str(seq).upper()
    # Get the standard codon table
    table = CodonTable.unambiguous_dna_by_id[1]

    # List all stop codons for this translation table
    stop_codons = table.stop_codons
    # print("Stop codons:", stop_codons)
    # Find positions of stop codons in the sequence
    stop_positions = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon in stop_codons:
            stop_positions.append(i)
    if len(stop_positions)>0:
        print("Warning: Stop codons found at positions:", stop_positions)


def validate_upstream_downstream(seq_str):
    # Convert to uppercase for consistency
    seq = Seq(seq_str.upper())
    
    # Condition 1: must be at least 50 bases long
    if len(seq) < 50:
        return False, f"Sequence too short (length {len(seq)} < 50)"
    
    # Condition 2: must contain only valid DNA nucleotides (A, T, G, C)
    valid_bases = {"A", "T", "G", "C"}
    if not set(seq).issubset(valid_bases):
        return False, "Sequence contains invalid characters (non-DNA bases)"
    
    return True, "Valid upstream/left overhang sequence"

# To check for High GC or High AT regions
def Validate_GC_AT(seq, min_length):
    """
    Find continuous high GC or high AT regions in a DNA sequence.
    
    Args:
        seq (str or Seq): DNA sequence
        min_length (int): minimum run length to report
    
    Returns:
        list of tuples: (region_type, start, end, sequence)
    """
    
    seq = str(seq).upper() # Ensure the sequence is a string and uppercase
    results = []
    
    # Regex for runs of G/C or A/T of length >= min_length
    gc_pattern = re.compile(rf"[GC]{{{min_length},}}")
    at_pattern = re.compile(rf"[AT]{{{min_length},}}")
    
    for match in gc_pattern.finditer(seq):
        results.append(("GC", match.start(), match.end(), match.group()))
    
    for match in at_pattern.finditer(seq):
        results.append(("AT", match.start(), match.end(), match.group()))

    
    return results

def Validate_Homoploymer_Stretches(seq, min_length):
     """
    Homopolymer stretches (> 5bp), like a polyA tail etc.
        Example: "AAAAA"
    Args:
        seq (str or Seq): DNA sequence
        min_length (int): minimum homopolymer length
    
    Returns:
        list of tuples: (base, start, end, stretch)
    """
     seq= str(seq).upper()
     results = []

     """ Homopolymer stretches are continuous sequences of the same DNA base, '
     'which are typically adenine (A), thymine (T), guanine (G), and cytosine (C). """
    
    #  a_pattern = r"A{min_length,}" # , means no upper limit
    #  t_pattern = r"T{min_length,}"
    #  G_pattern = r"G{min_length,}"
    #  C_pattern = r"C{min_length,}"

    # List of patterns for each base
     pattern_list = [ fr"A{{{min_length},}}",
        fr"T{{{min_length},}}",
        fr"G{{{min_length},}}",
        fr"C{{{min_length},}}"]

     pattern_names = ['A','T','G','C']

     for i, pattern in enumerate(pattern_list):
         compiled_pattern = re.compile(pattern)
         for match in compiled_pattern.finditer(seq):
             results.append((f"poly{pattern_names[i]} tail,",match.start(), match.end(), match.group()))

    
    #  print(results)
     return results

def Validate_Repeated_Fragments(seq,min_length=16):
    """
    Detects repeated DNA fragments longer than min_len bp.
    
    Args:
        seq (str or Seq): DNA sequence
        min_len (int): minimum repeat length (default 16)
    
    Returns:
        list of tuples: (fragment, positions)
            fragment: the repeated subsequence
            positions: list of start indices where it occurs
    """
    seq = str(seq).upper()  # Ensure the sequence is a string and uppercase
    repeats={}
    #sliding window approach for minimum length of min_length
    for j in range(min_length, len(seq)+1):
        for i in range(len(seq)-j+1):
            fragment=seq[i:i+j]
            if fragment in repeats:
                repeats[fragment].append(i)
            else:
                repeats[fragment]=[i]

    #get the repeated fragments
    repeated_fragments = [(frag, pos) for frag, pos in repeats.items() if len(pos) > 1]

    if len(repeated_fragments)>=1: print(repeated_fragments)

    return repeated_fragments


def Validate_GC_Content(seq, min_gc):
    """
    High GC Content (> 70%) means that in a DNA sequence, more than 70% of the bases are either G (guanine) or C (cytosine)

    Args:
        seq (str or Seq): DNA sequence
        min_gc (int): minimum GC content percentage to report for High GC content alert
    
    Returns:
        string: Alert message if GC content is high with the percentage value
    """

    seq=str(seq).upper()

    total_length_of_seq=len(seq)
    no_of_G=seq.count('G')
    no_of_C=seq.count('C')
    total_GC=no_of_G+no_of_C
    GC_content= (total_GC/total_length_of_seq)*100

    if(GC_content>min_gc):
        print(f'Warning: High GC content with {GC_content: .2f}%')
        return (f'Warning: High GC content with {GC_content: .2f}%')
    
    return None


def Validate_Tm_Values(seq,low_tm,high_tm,tm_method):
    """
    checks for very low or very high Tm values for the sequence.
    Args:
        seq (str or Seq): DNA sequence
        low_tm (int): minimum Tm value to report for Low Tm alert
        high_tm (int): maximum Tm value to report for High Tm alert
    
    Returns:
        string: Alert message of Tm value is too low or too high


    Package and Function: Biopython(https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html), provides 3 ways to calculate Tm.
    """
    seq=str(seq).upper()

    """Three ways to calculate Tm values using Biopython package."""
    # print('%0.2f' % mt.Tm_Wallace(seq))
    # print('%0.2f' % mt.Tm_GC(seq))

    # print('%0.2f' % mt.Tm_NN(seq))

    if(tm_method=='Wallace'):
        # Tm_Wallace is a simple formula based on the number of G/C and A/T pairs in the sequence.
        wallace=mt.Tm_Wallace(seq)
        if(wallace<low_tm):
            print(f'Warning: Low Wallace Tm value with {wallace: .2f}C')
            # return (f'Warning: Low Tm value with {wallace: .2f}C')
        elif(wallace>high_tm):
            print(f'Warning: High Wallace Tm value with {wallace: .2f}C')
            # return (f'Warning: High Tm value with {wallace: .2f}C')
        return 
    
    



    Tm =mt.Tm_NN(seq) #Tm_NN is implements the SantaLucia nearest-neighbor (NN) thermodynamic method.

    if(Tm<low_tm):
        print(f'Warning: Low Tm value with {Tm: .2f}C')
        return (f'Warning: Low Tm value with {Tm: .2f}C')
    elif(Tm>high_tm):
        print(f'Warning: High Tm value with {Tm: .2f}C')
        return (f'Warning: High Tm value with {Tm: .2f}C')
    
    return None


def Validate_Temperature_Difference(seq, temp_diff):
    """
    checks for The difference between forward and reverse primers should be within this temperature difference.
    Args:
        seq (str or Seq): DNA sequence
        temp_diff (int): maximum temperature difference between forward and reverse primers
    """

    reverse_seq=str(seq.reverse_complement()).upper() # Reverse complement of the sequence

    seq=str(seq).upper()

    Tm_forward=mt.Tm_NN(seq)

    
    Tm_reverse=mt.Tm_NN(reverse_seq)
    
    difference=abs(Tm_forward-Tm_reverse)



    if(difference>temp_diff):
        print(f'Warning: High Tm difference between forward and reverse primers with {difference: .2f}C')
        return (f'Warning: High Tm difference between forward and reverse primers with {difference: .2f}C')
    
    return None


def Validate_Hairpin_Formation(seq):
    """
    checks for Hairpin formation in the sequence.
    Args:
        seq (str or Seq): DNA sequence
    
    Returns:
        string: Alert message if hairpin structure is detected

    Note: Used primer3 package from the docs present at https://libnano.github.io/primer3-py/api/bindings.html
    """
    seq=str(seq).upper()
   
    #  only works for sequences <60bp
    if len(seq)>60:
        print('Sequence length exceeds 60bp for hairpin validation.')
        return None

    hairpin = primer3.bindings.calcHairpin(seq)

    if(hairpin.structure_found):
        print(f'Warning: Hairpin structure detected')
    
    
    return None


def get_codon_table_data(codon_table_file):
    data=pd.read_csv(codon_table_file, sep='\t',comment="#")
   #print(data.columns.tolist())

    # Sort by RSCU descending
    data_sorted = data.sort_values(by="RSCU", ascending=False)

    # Pick the first codon (highest RSCU) per amino acid
    best_codons = data_sorted.groupby("Amino acid").first().reset_index()    

    print("best codons are :" ,best_codons)

     

def Validations(sequences,tm_method):

    for seq in sequences:
        get_Stop_Codon(seq)
        min_length_gc_at = 8 # minimum length of continuous GC or AT to report
        GC_AT_result=Validate_GC_AT(seq, min_length_gc_at)
        if len(GC_AT_result)>=1: print(f"Warning: High GC or AT {GC_AT_result}")

        min_length_homopolymer = 6 # minimum homopolymer length
        
        Homoploymer_Stretches_result=Validate_Homoploymer_Stretches(seq, min_length_homopolymer)
        if len(Homoploymer_Stretches_result)>=1: print(Homoploymer_Stretches_result)

        Validate_Repeated_Fragments(seq,16)

        min_gc=70 # minimum GC content percentage to report for High GC content alert
        Validate_GC_Content(seq,min_gc)

        low_tm=60 # minimum Tm value to report for Low Tm alert
        high_tm=75 # maximum Tm value to report for High Tm alert
        Validate_Tm_Values(seq,low_tm,high_tm,tm_method) 

        temp_diff=5 # maximum temperature difference between forward and reverse primers
        Validate_Temperature_Difference(seq, temp_diff)

        Validate_Hairpin_Formation(seq)








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
    #parser.add_argument('-c', '--Codon_Table_File', default='Primer_Design/Codon_Table_Standard.csv', help="Codon Translation Table Mapping File")

    parser.add_argument('-c','--NCBI_Codon_Table_Value', default=1,help="NCBI Codon Translation Table ID Value") 
    # ORF_FILE = 'Primer_Design/HMT.txt' 
    parser.add_argument('-orf', '--ORF_File', default='Primer_Design/HMT.txt', help="File containing Open Reading Frame Sequence")

    # CODON_Statistics_FILE '
    parser.add_argument('-cod', '--Codon_Stat_File', default='nuclear_codon_statistics.tsv', help="File for Codon") 

     # Slecting the TM method '
    parser.add_argument('-tm', '--TM_Method', default='SantaLucia', help="Selecting the TM method") 

    args = parser.parse_args()

    out_path, path_fwd, path_rev = process_outputs(args.Output_Directory,args.Primer_Output_File,args.Forward_Primers_File,args.Reverse_Primers_File)
    

    remind_user_to_check_constants(args.Mutation_List,args.Output_Directory,args.ORF_File,args.NCBI_Codon_Table_Value)#args.Codon_Table_File)
    print(f'Working Directory: {os.getcwd()} \nProcessing...')
    
    orf_seq, mutations = read_orf_and_mutation_list(args.ORF_File,args.Mutation_List,args.NCBI_Codon_Table_Value)

    find_repeated_kmers(orf_seq)
    #check if provided mutations align with translation
    validate_mutations(mutations['Mutations'].tolist(), orf_seq) 
    print("Hi mutate ",mutations['Mutations'].tolist())

    primers = design_primers(orf_seq, mutations, args.NCBI_Codon_Table_Value)
    Validate_primer_length(primers["Length"].tolist())
    
    primer_order = create_primer_order_file(primers)
    primer_order.to_csv(out_path, index=False)
    

    sequences=primers["Sequence"].tolist() #list of primer sequences

    Validations(sequences, tm_method=args.TM_Method)
    
    selected_organism=args.Codon_Stat_File
    print(f"Selected organism for codon usage statistics: {selected_organism}")
    get_codon_table_data(f'Codon_Table/{selected_organism}')
    

    

    separate_primers_by_type(primer_order,path_fwd,path_rev)
    print(f"\nFinished in {time.time() - start_time:.6f} seconds.")