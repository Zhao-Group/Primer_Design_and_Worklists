import pandas as pd
import itertools
from collections import defaultdict, Counter
import sys, os
from pathlib import Path
import time


""" This script works for designing PCRs for Round2 .
"""

BASE_DIR = Path.cwd()
output = "Templates/Output_Phy_R2"
output_path = BASE_DIR / output
output_path.mkdir(parents=True, exist_ok=True)

Filename = "Phy_R2_processed.csv"
Out_file = output_path/Filename

Filename_2 = "Phy_Primers_and_R2mut.csv"
Out_file_2 = output_path/Filename_2
Filename_3 = "Echo_R2Phy_PCR.csv"
Out_file_3 = output_path/Filename_3
Filename_4 = "Phy_R2_Primers_Fluent.csv"
Out_file_4 = output_path/Filename_4

## PRIMER FILES for Round1 mutagenesis
Primer_plate1 = 'Templates/Phytase_Inputs/Phy_Plate1_Primers.csv'
Primer_plate2 = 'Templates/Phytase_Inputs/Phy_Plate2_Primers.csv'

Current_round_mutants = "Templates/Phytase_Round2.csv"


def Input_files():
    """Verify the input files"""

    """ Enter correct file names. """
    # For round 2, we need to find the most common mutant in Round 2
    # that will become the PCR template, and other mutants will use primer.
    
    print("\n⚠️  Reminder: Please check the following constants:")
    print(f"  - Current_mutants: {Current_round_mutants}")
    print(f"  - Output_file: {Out_file}")
    print("Modify these values at the top of the script if necessary.\n")
    
    df = Efficient_PCRdesign(Current_round_mutants)
    
    return df


def Efficient_PCRdesign(Current_round_mutants):
    
    Current_mutants = pd.read_csv(Current_round_mutants)
    
    Current_element_frequency = calculate_element_frequency(Current_mutants)
    df = process_frequency_data(Current_mutants,Current_element_frequency)
            
    ## This df will become input for next part of the script
    ## We will assign the PCR mutagenesis primers positions.
    ## Find out the location of tmeplates. 
    ## Create Fluent and Echo worklists. 
    
    return df


def calculate_element_frequency(Current_round_mutants):

    freq_list = Current_round_mutants["Mutation"].apply(
        lambda mut: list(itertools.permutations([m.strip() for m in mut.split('/')], 1))
    ).explode().tolist()  # Flatten the list of lists
    
    print(f"Total frequency entries: {len(freq_list)}")
    
    return Counter(freq_list)



def Get_primer_wells():
    df1 = pd.read_csv(Primer_plate1)
    df2 = pd.read_csv(Primer_plate2)

    # Initialize a defaultdict to store mutations and their corresponding wells
    mutation_dict = defaultdict(list)
    def update_dict(df):
        for _, row in df.iterrows():
            mutation = row['Mutation']
            well = str(row['Well'])
            mutation_dict[mutation].append(well)
    
    update_dict(df1)
    update_dict(df2)
    
    return dict(mutation_dict)


def process_frequency_data(Current_mutants, Current_element_frequency):
    """Process Quad data to find matching triple mutations."""
    Curr_list, Curr_well = [], [] # these will contain the mutants from previous round
    
    well_dict = Get_primer_wells()
    for _, row in Current_mutants.iterrows():
        split_mutants = list(map(str.strip, row["Mutation"].split('/')))
        mut_combs = itertools.permutations(split_mutants, 1)
        ordered_combs = sorted(mut_combs, key=lambda x: Current_element_frequency[x], reverse=True)
        item = next((i for i in ordered_combs if i in Current_element_frequency), None)

        select_mutation = list(item)[0]
        
        try:
            if select_mutation:
                Curr_well.append(well_dict[select_mutation][0])
                Curr_list.append(select_mutation)
            else:
                print("WARNING! For some mutations, no templates were found.")
                Curr_list.append('None')    
                Curr_well.append('None')
                
        except KeyError:
                Curr_list.append('Wild type')    
                Curr_well.append('TBD')

    
    Current_mutants["PCR_template"] = Curr_list
    Current_mutants["PCR_template_well"] = Curr_well

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
    print(f'\nWorking Directory: {os.getcwd()} \nProcessing...')
    
    df = Input_files()
    df = Start_Worklist_Scripts(df)
    df.to_csv(Out_file, index=False)

    print(f"\nFinished in {time.time() - start_time:.6f} seconds.")
