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



def find_optimal_codon(orf_seq, pos, codons):
    """Find the most similar codon to the original one."""
    original_codon = str(orf_seq[pos * 3:(pos + 1) * 3])
    return max(codons, key=lambda codon: sum(a == b for a, b in zip(codon, original_codon)))

def validate_position(amino_acid_pos, protein_len):
    """Ensure amino acid position is valid and within the ORF."""
    if not (1 <= amino_acid_pos <= protein_len):
        sys.exit('Amino acid position exceeds ORF bounds.')

def extract_primer(seq, start):
    """Extract a 27bp primer around the mutation site."""
    return seq[max(0, start - 12):min(len(seq), start + 15)]



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


def find_repeated_kmers(seq, k=16):
    """Check for repeated k-mers in the sequence. Change the k-value to check for 
    smaller repeats if needed."""
    seq_str = str(seq)  # Ensure we're working with a string
    kmers = [seq_str[i:i+k] for i in range(len(seq_str) - k + 1)]
    repeats = sum(seq_str.count(kmer) > 1 for kmer in kmers)
    
    if repeats:
        print(f'\nWarning: {repeats} repeats of {k}bp or more in the ORF. This will affect PCR & SDM.')