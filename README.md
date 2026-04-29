RMSDistance

This repository contains the analysis scripts used for the publication “Ligand-induced activation of RyR1 in native membranes.” The scripts reproduce the key computational steps used in the study and are provided for transparency and reproducibility.

## Contents

- **ChimeraX .py scripts** for running structure-related analysis in no-GUI mode.
- **MATLAB .m script** for the multidimensional scaling.

## Requirements
- Python 3.8.10 or later
- UCSF ChimeraX 1.8 or later
- MATLAB R2021a or later
- Any input files referenced by the scripts, like the model list .txt file which is provided as a reference.

## Usage

### ChimeraX
Run the ChimeraX scripts in command-line / no-GUI mode using your local ChimeraX installation.

1. fetch_changechains_align_save.py
this will fetch all the PDB models from the database, align it to the reference model and update the chain IDs to match the reference model (A,B,C,D...);
all the modified models will be saved to the ./models

run cmd example:
chimerax --nogui --script fetch_changechains_align_save.py > fetch_changechains_align_save.log 2>&1

3. universal_align_rmsd_normalized_all.py
this will calculate the normalized RMSD pairwise for all of the provided models and for all of the provided reference and target domain combinations.
### MATLAB
Open and run the MATLAB scripts in MATLAB with the required input data in place.
3. rmsd.m
this will take the output files from the pairwise RMSD calculation done previously, assemble it into a single distance matrix and perform a non-metric multidimensional scaling on it.
the results can be used to create a 2D plot.

## Notes

- File paths may need to be updated to match your system.
- Some scripts may depend on external data or model files not included here.
- These scripts were used to generate results reported in the associated publication.
