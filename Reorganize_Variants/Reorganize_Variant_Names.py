import pandas as pd
from collections import Counter
from pathlib import Path

# re organize the mutants in a way that each mutant contains 
# most common mutant first and so on...

def Re_organize_mutation_names(Output_DIR, File_dir):
    # Get all .csv files in the subdirectory
    input_files = list(File_dir.glob("*.csv"))

    # Iterate over all matching files
    for file in input_files:
        # Load the data from each file
        quad_df = pd.read_csv(file)
        # Create a flat list of all mutations across rows
        all_mutants = [
            mut.strip() 
            for mutation in quad_df["Mutation"].astype(str) 
            for mut in mutation.split('/')
        ]

        # Compute frequency of each mutation
        mutation_frequency = Counter(all_mutants)
        # Reorder mutations in each row based on frequency
        quad_df["Mutation"] = quad_df["Mutation"].apply(
            lambda x: ' / '.join(
                sorted(x.split('/'), key=lambda m: mutation_frequency[m.strip()], reverse=True)
            )
        )

        new_filename = file.stem + "_reorganized.csv"
        new_filepath = Output_DIR / new_filename
        quad_df.to_csv(new_filepath, index=False)
        print(f"Processed and saved: {new_filename}")
        
if __name__ == "__main__":

    Output_DIR = Path.cwd() / "Output_Reorganized"
    Output_DIR.mkdir(parents=True, exist_ok=True)

    # Define the subdirectory that conatins all files 
    File_dir =  Path.cwd()

    Re_organize_mutation_names(Output_DIR, File_dir)