import pandas as pd
import itertools
from collections import defaultdict, Counter
import sys, os
from pathlib import Path
import time


""" 
This script works for designing PCRs for Round3 onwards.
For the current mutants in Round 4, the previous mutants file will be Round 3 one. 
For Round 2 PCRs, there is a separate script. 
"""

BASE_DIR = Path.cwd()
output = "Templates/Output_HMT_R3"
output_path = BASE_DIR / output
output_path.mkdir(parents=True, exist_ok=True)


Filename = "HMT_R3_processed.csv"
Out_file = output_path/Filename
Filename_2 = "Primers_and_R3mut_forR4.csv"
Out_file_2 = output_path/Filename_2
Filename_3 = "Echo_R4HMT_PCR.csv"
Out_file_3 = output_path/Filename_3
Filename_4 = "HMT_R4_Primers_Fluent.csv"
Out_file_4 = output_path/Filename_4


k = 2       ## For 3 mutations k = 2, for 4 mutations k = 3 and so on

## PRIMER FILES for generating worklists for SDM
Primer_plate1 = 'Templates/HMT_Inputs/HMT_Primers_Plate1.csv'
Primer_plate2 = 'Templates/HMT_Inputs/HMT_Primers_Plate2.csv'

# the organized file with mutants/well from previous round
Prev_round_mutants = "Templates/HMT_Round2.csv"  
# the organized file with mutants from current round, for which we need 
# we need PCR templates from previous round.
Current_round_mutants = "Templates/HMT_Round3.csv"


def Input_files():
    """
    Verify the input files
    """
    print("\n⚠️  Reminder: Please verify PATHS and other CONSTANTS:")
    print(f"  - Previous_mutants: {Prev_round_mutants}")
    print(f"  - Current_mutants: {Current_round_mutants}")
    print(f"  - Output_file: {Out_file}")
    print("There are more at the top of the script.\n")
    
    df = Efficient_PCRdesign(Prev_round_mutants, Current_round_mutants)
    
    return df


def Efficient_PCRdesign(Prev_round_mutants,Current_round_mutants):
    
    Prev_mutants = pd.read_csv(Prev_round_mutants)
    # expand the n object to tuple by organizing those by their frequency
    Prev_list, well_dict = prepare_prev_data(Prev_mutants)
    
    Current_mutants = pd.read_csv(Current_round_mutants)
    # all possible combination of 3 object for 4 object, as a list
    # Maximum allowable is 8 object with 7 permutations 
    
    Current_element_frequency = calculate_element_frequency(Current_mutants, k)
    df = process_frequency_data(Current_mutants, Prev_list, well_dict, Current_element_frequency,k)
    
    return df


def prepare_prev_data(Prev_mutants):
    """Prepare R3 mutants data."""
    Prev_list, well_dict = [], {}
    for _, row in Prev_mutants.iterrows():
        split_mutants = tuple(map(str.strip, row["Mutation"].split('/')))
        Prev_list.append(split_mutants)
        well_dict[split_mutants] = row["Well"]
    return Prev_list, well_dict


def calculate_element_frequency(Current_round_mutants, k):
    """Calculate frequency of mutant combinations."""
    Mutation_lengths = Current_round_mutants["Mutation"].apply(lambda x: len([m.strip() for m in x.split('/')]))
    if (max(Mutation_lengths) > 8 or len(set(Mutation_lengths)) > 1):
        sys.exit("Either more than 8 mutations or Incosistent length!!")

    freq_list = Current_round_mutants["Mutation"].apply(
        lambda mut: list(itertools.permutations([m.strip() for m in mut.split('/')], k))
    ).explode().tolist()  # Flatten the list of lists
    
    print(f"Total frequency entries: {len(freq_list)}")
    
    return Counter(freq_list)


def process_frequency_data(Current_mutants, Prev_list, well_dict, Current_element_frequency, k):
    """Process Quad data to find matching triple mutations."""
    Curr_list, Curr_well = [], [] # these will contain the mutants from previous round

    for _, row in Current_mutants.iterrows():
        split_mutants = list(map(str.strip, row["Mutation"].split('/')))
        mut_combs = itertools.permutations(split_mutants, k)
        ordered_combs = sorted(mut_combs, key=lambda x: Current_element_frequency[x], reverse=True)

        item = next((i for i in ordered_combs if i in Prev_list), None)
        
        if item:
            Curr_list.append(' / '.join(item))
            Curr_well.append(well_dict[item])
        else:
            print("WARNING! For some mutations, no templates were found.")
            Curr_list.append('None')    
            Curr_well.append('None')
            
    Current_mutants["PCR_template"] = Curr_list
    Current_mutants["PCR_template_original_well"] = Curr_well

    return Current_mutants


#############################################################
""" Second part of the script for worklist generation. """
#############################################################


def Start_Worklist_Scripts(df):
    
    quad_df = df

    single_mutant_list = extract_single_mutants(df)
    df["Single"] = single_mutant_list
    
    well96 = well96_list()
    df["PCR_Plate"] = well96[:len(df)]
    
    single_mut_well = add_single_mutant_well(df)
    df["Primer Well"] = single_mut_well
        
    #df.to_csv(Out_file, index=False) # Save this file as main processed file for reference
    
    df = prepare_for_next_round(df)
    
    return df
    

def extract_single_mutants(Processed_df):

    """Extract single mutants, which needs mutagenesis primers."""
    single_mut_list = [
        mutant.strip()
        for _, row in Processed_df.iterrows()
        for mutant in row["Mutation"].split('/')
        if mutant.strip() not in map(str.strip, row["PCR_template"].split('/'))
    ]
    
    return single_mut_list


def well96_list():
    """Define PCR plate wells."""
    rows, cols = 'ABCDEFGH', range(1, 13)
    well96 = [f"{row}{col}" for row in rows for col in cols]
    
    return well96


def add_single_mutant_well(quad_df):
    """Add the location of mutagenesis primers."""
    primer_wells = {}
    
    ## Update the location/name of primer files
    primer_wells.update(load_plate_data(Primer_plate1, "P1"))
    primer_wells.update(load_plate_data(Primer_plate2, "P2"))

    single_mut_well = [
        primer_wells.get(str(row["Single"]).strip(), '') for _, row in quad_df.iterrows()
    ]
    
    return single_mut_well


def load_plate_data(filename, prefix):
    """Load primer plate data and add prefix."""
    plate_data = pd.read_csv(filename)
    return {str(row["Mutation"]).strip(): f"{prefix} {row['Well']}" for _, row in plate_data.iterrows()}


def prepare_for_next_round(df):
    """Prepare data for the next round."""
    quad_df = df
    frequency_counts = df["PCR_template"].value_counts()

    new_df = pd.DataFrame({
        "PCR_template": frequency_counts.index,
        "R3 Well": [ 
            df.loc[df["PCR_template"] == PCR_template, "Well"].iloc[0]
            for PCR_template in frequency_counts.index
        ],
        "Count_I": frequency_counts.values,
        "RE_96Well": well96_list()[:len(frequency_counts)],
        "Template_384Well": assign_384well_positions(len(frequency_counts)),
        "Blank": ''
    })

    single_counts = df["Single"].value_counts()
    single_df = pd.DataFrame({
        "Single": single_counts.index,
        "Primer Well": [
            quad_df.loc[df["Single"] == single, "Primer Well"].iloc[0]
            for single in single_counts.index
        ],
        "Count_II": single_counts.values
    })

    df = assign_template_positions_to_quad(df, new_df)
    new_df = pd.concat([new_df, single_df], axis=1, ignore_index=False)
    
    new_df.to_csv(Out_file_2, index=False)

    make_worklists(df, new_df)
    
    return df
    
        
def assign_384well_positions(length):
    """Assign positions in a 384-well plate."""
    rows, cols = 'ABCD', range(1, 25)
    return [f"{row}{col}" for row in rows for col in cols][:length]


def assign_template_positions_to_quad(df, new_df):
    
    """Assign template positions to Quad data in dictionary."""
    template_384well = {
        row["PCR_template"]: row["Template_384Well"]
        for _, row in new_df.iterrows()
    }
    
    df["PCR_template_Echo_384Well"] = [
        template_384well.get(PCR_template, '') for PCR_template in df["PCR_template"]
    ]
    
    columns = df.columns[:-1].tolist()  
    columns.insert(4, "PCR_template_Echo_384Well") 
    df_rearranged = df[columns]
    
    return df_rearranged


def make_worklists(df, new_df):
    
    """Create worklists for Echo and Fluent systems."""
    echo_df = df[["PCR_template_Echo_384Well", "PCR_Plate", "Mutation"]].copy()
    echo_df["Volume"] = 2500    
    echo_df = echo_df[["PCR_template_Echo_384Well", "PCR_Plate", "Volume", "Mutation"]]
    echo_df.to_csv(Out_file_3, index=False)

    fluent_df = pd.DataFrame({
        "Source Plate": df["Primer Well"].str.split().str[0].replace({"P1": "Plate1", "P2": "Plate2"}),
        "Source Well": df["Primer Well"].str.split().str[1],
        "Dest Plate": "PCR_Plate",
        "Dest Well": df["PCR_Plate"],
        "Volume": 12.5
    })
    fluent_df.to_csv(Out_file_4, index=False)


if __name__ == "__main__":
    
    start_time = time.time()
    print(f'Working Directory: {os.getcwd()} \nProcessing...')

    df = Input_files()
    df = Start_Worklist_Scripts(df)
    df.to_csv(Out_file, index=False)    
    
    print(f"\nFinished in {time.time() - start_time:.6f} seconds.")