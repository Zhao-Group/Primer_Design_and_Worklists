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
    """Extract a primer window around the mutation site.

    The -12/+15 window (12bp upstream of the codon, codon + 12bp downstream)
    is the method's fixed geometry. When overhangs flank the ORF the window can
    extend into them, giving full-length primers for terminal mutations.
    TODO(SME q2, 2026-06-22): confirm the flank size / window geometry; see
    MUTAGENESIS_INTEGRATION_SPEC.md.
    """
    return seq[max(0, start - 12):min(len(seq), start + 15)]



def design_primers(orf_seq, mutations, codon_table_id, left_overhang='', right_overhang=''):
    """Main primer design function.

    Primers are windowed over `left_overhang + ORF + right_overhang` so that
    mutations near the ORF termini still yield full-length primers when overhangs
    are supplied. Mutation positions are validated/translated against the ORF
    alone; only the windowing is offset by len(left_overhang).
    """
    primers = []
    left = str(left_overhang or '').upper()
    right = str(right_overhang or '').upper()
    offset = len(left)
    full = Seq(left) + orf_seq + Seq(right)
    protein_len = len(orf_seq) // 3

    for mutation in mutations['Mutations']:
        mutation = str(mutation).strip()
        pos = int(re.search(r'\d+', mutation).group()) - 1
        validate_position(pos + 1, protein_len)

        # Get all possible codons coding for the target aa residue
        original_aa, target_aa = mutation[0], mutation[-1]
        codon_table = CodonTable.unambiguous_dna_by_id[int(codon_table_id)].forward_table
        codons = [key for key, value in codon_table.items() if value == target_aa]
        if not codons:
            sys.exit(f'ERROR: target residue "{target_aa}" in mutation "{mutation}" is '
                     f'not a valid amino acid for codon table {codon_table_id}.')
        new_codon = find_optimal_codon(orf_seq, pos, codons)

        codon_start = offset + pos * 3
        mutated_seq = full[:codon_start] + new_codon + full[codon_start + 3:]
        primer = extract_primer(mutated_seq, codon_start)
        tm = int(math.ceil(calc_tm(str(primer), dv_conc=2, tm_method='santalucia', salt_corrections_method='owczarzy')))

        primers.append((mutation, str(primer), tm, int(SeqUtils.gc_fraction(primer)*100.0), len(primer), pos + 1)) #SeqUtils.GC is deprecated after Biopython 1.82

    return pd.DataFrame(primers, columns=['Name', 'Sequence', 'Tm', 'GC', 'Length', 'AA_Position'])


def find_repeated_kmers(seq, k=16):
    """Check for repeated k-mers in the sequence. Change the k-value to check for 
    smaller repeats if needed."""
    seq_str = str(seq)  # Ensure we're working with a string
    kmers = [seq_str[i:i+k] for i in range(len(seq_str) - k + 1)]
    repeats = sum(seq_str.count(kmer) > 1 for kmer in kmers)
    
    if repeats:
        print(f'\nWarning: {repeats} repeats of {k}bp or more in the ORF. This will affect PCR & SDM.')