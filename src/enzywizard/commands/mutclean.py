from __future__ import annotations
from argparse import Namespace
from ..services.mutclean_service import run_mutclean_service

def add_mutclean_parser(subparsers) -> None:
    parser = subparsers.add_parser("mutclean",help="Clean a pair of wild-type and mutant CIF/PDB structure files with a specified amino acid substitution.")
    parser.add_argument("-w","--wt_input_path", required=True, help="Path to wild-type protein structure file (CIF/PDB).")
    parser.add_argument("-m","--mut_input_path", required=True, help="Path to mutant protein structure file (CIF/PDB).")
    parser.add_argument("-s","--mutation", required=True, help="Amino acid substitution(s) describing mutations, e.g., A123V or A123V/G456D.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to a directory for outputting cleaned CIF, PDB, and FASTA files and a JSON report.")
    parser.add_argument( "--add_H", type=lambda x: str(x).lower() in ["true", "1", "yes"], default=True, help="Whether to add hydrogens using OpenMM (True/False, default: True)." )
    parser.add_argument("--pH",type=float,default=7.0,help="pH value for hydrogen addition (default: 7.0).")
    parser.set_defaults(func=run_mutclean)

def run_mutclean(args: Namespace) -> None:
    run_mutclean_service(wt_input_path=args.wt_input_path, mut_input_path=args.mut_input_path,mutation=args.mutation,output_dir=args.output_dir,add_H=args.add_H, pH=args.pH)

# input parameters:
'''
-w --wt_input_path Required. Path to the wild-type protein structure file (CIF or PDB).
-m --mut_input_path Required. Path to the mutant protein structure file (CIF or PDB).
-s --mutation Required. Amino acid substitution(s), e.g.: A123V A123V/G456D Multiple mutations should be separated by '/'.
-o --output_dir Required. Output directory for cleaned structures, sequences, and report.
--add_H Optional. Whether to add hydrogens using OpenMM (default: True).
--pH Optional. pH value used for hydrogen addition (default: 7.0).
'''

# output content:
'''
The program outputs:

1. Cleaned CIF and PDB files for both wild-type and mutant structures:
    - cleaned_<wt_name>.cif
    - cleaned_<wt_name>.pdb
    - cleaned_<wt_name>.fasta
    - cleaned_<mut_name>.cif
    - cleaned_<mut_name>.pdb
    - cleaned_<mut_name>.fasta
2. A JSON report containing:

    - "output_type": "enzywizard_mutclean"

    - "amino_acid_substitution":
        The original input mutation(s).

    - "cleaned_amino_acid_substitution":
        Mutation(s) remapped to cleaned residue indices.

    - "wt_amino_acid_mapping_old_to_new":
        Mapping from original residues to cleaned residues for wild-type.
        Each entry includes:
            - original residue index and name
            - cleaned residue index and name
            - hydrogen atom counts (before and after cleaning)

    - "wt_clean_statistics":
        Cleaning statistics for wild-type, including:
            - residue name standardization count
            - removal of non-standard residues
            - removal of residues missing backbone atoms
            - removal of residues missing required heavy atoms
            - removal of residues with invalid occupancy
            - removal of insertion codes
            - total number of retained residues

    - "mut_amino_acid_mapping_old_to_new":
        Same mapping information for mutant structure.

    - "mut_clean_statistics":
        Cleaning statistics for mutant structure.
'''

# functionality/process:
'''
The program processes a pair of wild-type and mutant protein structures through the following steps:

1. Input validation:
    - Check that both input files exist
    - Ensure filenames are valid and distinct
    - Validate mutation format and positions against sequence lengths

2. Structure loading:
    - Load both structures using Biopython
    - Extract a single chain from each structure

3. Structure cleaning (performed independently for wild-type and mutant):
    - Keep only protein residues (remove hetero residues)
    - Standardize residue names using MODRES mapping
    - Remove non-standard residues
    - Remove residues missing backbone atoms (N, CA, C)
    - Remove residues missing required heavy atoms
    - Remove residues with invalid occupancy
    - Remove insertion codes
    - Renumber residues consecutively starting from 1
    - Build a mapping from original residues to cleaned residues

4. Mutation remapping:
    - Map original mutation positions to cleaned residue indices
    - Ensure consistency between wild-type and mutant mappings

5. Clean structure validation (performed before hydrogen addition):
    - Ensure the cleaned structure:
        - contains exactly one model and one chain (chain ID = 'A')
        - has continuous residue numbering
        - contains only standard amino acids
        - has complete backbone and required heavy atoms
    - Verify that residue coordinates (e.g., CA atoms) are preserved after cleaning

6. Optional hydrogen addition:
    - Add hydrogens using OpenMM Modeller with specified pH and force field
    - Convert the resulting structure back to Biopython format

7. Output generation:
    - Save cleaned structures in CIF and PDB formats and sequences in FASTA format
    - Generate a JSON report summarizing mappings and cleaning statistics
'''

# dependency:
'''
Biopython
OpenMM (for optional hydrogen addition)
MODRES residue mapping from:
https://www.wwpdb.org/data/ccd
'''
# reference:
'''
https://docs.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-structures
'''


