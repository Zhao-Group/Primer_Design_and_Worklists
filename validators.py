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

def Validate_primer_length(primer_length):
    
    count = sum(1 for x in primer_length if x < 27)
    # If the count is greater than 0, print "Yes"
    if count:
        print(f'WARNING: {count} primers less than 27bp in the file.')

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