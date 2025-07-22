## Primer design, PCR organization, and worklist generation
![Summary]( Figures/Schema.png)

### The Manuscript
This repository accompanies the work [A generalized platform for artificial intelligence-powered autonomous enzyme engineering
](https://www.nature.com/articles/s41467-025-61209-y)
<br>Please cite the reference below if you use this tool.

###  PCR Primer Design and Variant Management
This repository contains a set of Python scripts designed for efficient primer design and PCR template selection across multiple rounds of continuous machine-learning guided autonomous protein engineering. Each script handles specific tasks related to mutagenesis primer design, variant naming, and PCR preparations for successive rounds.

### Repository Contents
- **Primer_Design_SDM.py**: Script for designing primers for site-directed mutagenesis (SDM) using variant names. This script needs a text file with ORF sequence and a CSV file with name of mutations. The input files need to be divided as 96 mutations per file for 96-well throughput assay. The output can be used to directly order primers in 96-well plate. <br><br>
- **Reorganize_Variant_Names.py**: This script reorganizes variant names based on frequency. It reads a CSV file of mutations predicted by the supervised learning models and sorts the mutation names according to their occurrence. <br><br>
- **Round2_PCR_Templates_Selection.py**: This script is specific for second round of PCR. For the PCR, the variants from previous round are used as PCR templates for next round. This script minimizes the number of templates needed, ensuring efficient PCR design. It processes input files to determine the most common variants, assigns wells in 96-well PCR plate, and generates worklists for PCR for Echo and Fluent systems. <br><br>
- **Round3_Onwards_PCRTemplate_Selection.py**: This script is similar to above but is tailored for round 3 and subsequent rounds. It calculates template requirements based on prior rounds' mutations, processes template selections, and generates worklists. 

###  Important Notes
- Each script has example input files and output files in associated folders. 
- Verify constants at the beginning of the script for the correct file paths.
- Each script expects certain directories and files to be present.

###  Requirements
- Python 3.x
- Libraries: `pandas`, `itertools`, `collections`, `pathlib`, `os`, `sys`, `time`

### Reference
<details>
<summary>If you use this tool, please cite us:</summary>

```bibtex
Singh, Nilmani, et al. "A generalized platform for artificial intelligence-powered autonomous enzyme engineering" Nature Communications, vol. 16:5648 (2025).
```
</details>
