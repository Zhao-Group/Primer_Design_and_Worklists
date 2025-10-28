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

def get_codon_table_data(codon_table_file):
    data=pd.read_csv(codon_table_file, sep='\t',comment="#")
   #print(data.columns.tolist())

    # Sort by RSCU descending
    data_sorted = data.sort_values(by="RSCU", ascending=False)

    # Pick the first codon (highest RSCU) per amino acid
    best_codons = data_sorted.groupby("Amino acid").first().reset_index()    

    print("best codons are :" ,best_codons)