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

# Create the full path to the file ###########################
######################################
def process_outputs(output_dir,primer_output_file,fwd_primers_file,rev_primers_file):
    base_dir=Path.cwd()
    Output_Path = base_dir / output_dir / primer_output_file
    Output_Path.parent.mkdir(parents=True, exist_ok=True)
    Output_Path_Fwd = base_dir / output_dir / fwd_primers_file
    Output_Path_Rev = base_dir / output_dir / rev_primers_file
    return Output_Path, Output_Path_Fwd, Output_Path_Rev

def parse_mutation(mutation_str):
    """
    Parses a mutation string into type, position, original, and target amino acids.
    
    Supports:
        Substitution: M78Y
        Deletion:     M78del
        Insertion:    M78insY
    Returns:
        dict: {type, position (1-based), original_aa, new_aa}
    """
    mutation_str = mutation_str.strip()

    # Substitution: e.g., M78Y
    sub_match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", mutation_str)
    if sub_match:
        return {
            "type": "substitution",
            "position": int(sub_match.group(2)),
            "original_aa": sub_match.group(1),
            "new_aa": sub_match.group(3)
        }

    # Deletion: e.g., M78del
    del_match = re.fullmatch(r"([A-Z])(\d+)del", mutation_str)
    if del_match:
        return {
            "type": "deletion",
            "position": int(del_match.group(2)),
            "original_aa": del_match.group(1),
            "new_aa": None
        }

    # Insertion: e.g., M78insY
    ins_match = re.fullmatch(r"([A-Z])(\d+)ins([A-Z]+)", mutation_str)
    if ins_match:
        return {
            "type": "insertion",
            "position": int(ins_match.group(2)),
            "original_aa": ins_match.group(1),
            "new_aa": ins_match.group(3)
        }

    raise print(f"Invalid mutation format: {mutation_str}")


def read_orf_and_mutation_list(orf_file,mutation_list_file,codon_table_file_id):
    """Read the ORF sequence, mutation list, and codon table."""
    orf_seq = Seq(Path(orf_file).read_text().strip())

    # ✅ Check if cDNA starts with start codon ATG
    if not orf_seq.upper().startswith("ATG"):
        print("⚠️ Warning: The provided cDNA sequence does not start with ATG (start codon).")
        
    mutation_list = pd.read_csv(mutation_list_file)

     # Limit check for mutation rows
    if len(mutation_list) > 1000:
        print(f"❌ ERROR: Mutation list contains {len(mutation_list)} entries — limit is 1000 per run. Please split the file and retry.")
    elif len(mutation_list) == 0:
        print("❌ ERROR: Mutation list is empty. Please provide at least one mutation.")
    #codon_table = CodonTable.unambiguous_dna_by_id[codon_table_file_id].tolist() #pd.read_csv(codon_table_file)
    return orf_seq, mutation_list #, codon_table


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