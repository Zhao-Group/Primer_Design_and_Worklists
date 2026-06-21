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
import json
from pathlib import Path
import argparse

from validators import get_primer_alerts

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


def read_overhang(path):
    """Read an optional overhang sequence file; '' if not provided/missing."""
    if not path:
        return ''
    p = Path(path)
    if not p.exists():
        return ''
    return p.read_text().strip()


# --- Plate layout -----------------------------------------------------------
PLATE_ROWS = "ABCDEFGH"
PLATE_COLS = 12
WELLS_PER_PLATE = len(PLATE_ROWS) * PLATE_COLS  # 96

# TODO(SME q1, 2026-06-22): plate topology unconfirmed. Defaults below; flip the
# two constants once the SME confirms. See MUTAGENESIS_INTEGRATION_SPEC.md.
#  - forward & reverse primers go on SEPARATE plate sets
#  - within a plate, fill row-major (A1..A12, B1..)
#  - a forward/reverse pair shares the same well index on its respective plate
#  - SME-CONFIRMED: ordered by aa position, variants at a site grouped together
FORWARD_REVERSE_SEPARATE_PLATES = True
FILL_ORDER = "row-major"  # or "column-major"

# Output columns for the single combined results CSV consumed by the backend.
RESULT_COLUMNS = ['Mutations', 'Assembly_Fragments', 'Direction', 'Sequence',
                  'Tm', 'GC', 'Length', 'Alerts', 'Plate', 'Well']


def _plate_wells():
    if FILL_ORDER == "column-major":
        return [f"{r}{c}" for c in range(1, PLATE_COLS + 1) for r in PLATE_ROWS]
    return [f"{r}{c}" for r in PLATE_ROWS for c in range(1, PLATE_COLS + 1)]


def _make_row(src, direction, seq, plate, well, tm_method):
    """One combined-CSV row. Alerts are JSON-encoded (a list) so the backend can
    json.loads them straight into the typed result."""
    return {
        'Mutations': src['Name'],
        'Assembly_Fragments': 1,  # v1 single-substitution: one fragment per mutation
        'Direction': direction,
        'Sequence': seq,
        'Tm': src['Tm'],
        'GC': src['GC'],
        'Length': len(seq),
        'Alerts': json.dumps(get_primer_alerts(seq, tm_method)),
        'Plate': plate,
        'Well': well,
    }


def build_results(primers, tm_method='SantaLucia'):
    """Build the single combined results table.

    One forward + one reverse row per designed primer, ordered by aa position
    (variants at the same site grouped), assigned across 96-well plates, each
    row carrying its own per-primer alerts.
    """
    wells = _plate_wells()
    ordered = primers.sort_values('AA_Position', kind='stable').reset_index(drop=True)
    n = len(ordered)
    fwd_plate_count = math.ceil(n / WELLS_PER_PLATE) if n else 0

    rows = []
    for i, src in ordered.iterrows():
        well = wells[i % WELLS_PER_PLATE]
        plate_offset = i // WELLS_PER_PLATE
        fwd_seq = src['Sequence']
        rev_seq = str(Seq(fwd_seq).reverse_complement())

        rows.append(_make_row(src, 'Forward', fwd_seq,
                              f"Plate {plate_offset + 1}", well, tm_method))

        rev_plate = (fwd_plate_count + plate_offset + 1
                     if FORWARD_REVERSE_SEPARATE_PLATES else plate_offset + 1)
        rows.append(_make_row(src, 'Reverse', rev_seq,
                              f"Plate {rev_plate}", well, tm_method))

    print(f"\nBuilt {len(rows)} primers across {fwd_plate_count} forward plate(s).")
    return pd.DataFrame(rows, columns=RESULT_COLUMNS)