from pathlib import Path
from Bio import SeqUtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from primer3 import calc_tm # as calcTm - calcTm is deprecated
import primer3

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Data import IUPACData

import pandas as pd
from itertools import combinations
import time
import os
import sys
import re
import math
from pathlib import Path
import argparse

from parsers import process_outputs, read_orf_and_mutation_list, create_primer_order_file, separate_primers_by_type
from design import design_primers, find_repeated_kmers
from validators import validate_mutations, Validate_primer_length, Validations
from codon_utils import get_codon_table_data



def remind_user_to_check_constants(mut_file,out_path,orf_file,codon_file_id):
    """Reminds the user to review constants and make changes if needed."""
    print("\n⚠️ Reminder: Please check the following inputs for correctness: \n")
    print(f"  - MUTATION_FILE: {mut_file}\n")
    print(f"  - PRIMER_OUTPUT: {out_path}\n")
    print(f"  - ORF_FILE: {orf_file}\n")
    print(f"  - CODON_TABLE: {CodonTable.unambiguous_dna_by_id[int(codon_file_id)]}\n")
    print("Modify these and other relevant values in CLI arguments if necessary. Run with -h for help.\n")


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