from __future__ import annotations
from argparse import Namespace
from ..services.mut_clean_service import run_mutclean_service

def add_mutclean_parser(subparsers) -> None:
    parser = subparsers.add_parser("mut_clean",help="Clean both a wild-type protein structure and its corresponding mutant protein structure, trace cleaned amino acid substitution(s), generate multi-format cleaned protein files (CIF, PDB, and FASTA), and provide a detailed traceable cleaning report.")
    parser.add_argument("-w","--wt_input_path", required=True, help="Path to wild-type protein structure file (CIF/PDB).")
    parser.add_argument("-m","--mut_input_path", required=True, help="Path to mutant protein structure file (CIF/PDB).")
    parser.add_argument("-s","--mutation", required=True, help="Amino acid substitution(s) describing mutations, e.g., A123V or A123V,G456D.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to a directory for outputting cleaned CIF, PDB, and FASTA files and a JSON report.")
    parser.add_argument("--no_add_H",action="store_false",dest="add_H",help="Disable adding hydrogens using OpenMM (default: enabled).")
    parser.set_defaults(add_H=True)
    parser.add_argument("--pH",type=float,default=7.0,help="pH value for hydrogen addition (default: 7.0).")
    parser.set_defaults(func=run_mutclean)

def run_mutclean(args: Namespace) -> None:
    run_mutclean_service(wt_input_path=args.wt_input_path, mut_input_path=args.mut_input_path,mutation=args.mutation,output_dir=args.output_dir,add_H=args.add_H, pH=args.pH)

# ==============================
# Command: enzywizard-mut-clean
# ==============================

# brief introduction:
'''
EnzyWizard-Mut-Clean is a command-line tool for cleaning both a wild-type
protein structure and its corresponding mutant protein structure, tracing 
cleaned amino acid substitution(s), generating
multi-format cleaned protein files (CIF, PDB, and FASTA), and providing a detailed 
traceable cleaning report. It standardizes residue
names, removes problematic residues (non-standard residues, residues with
missing backbone atoms, residues with missing required heavy atoms, residues
with unexpected heavy atoms, residues with invalid occupancy), repairs residue
order by renumbering residues continuously, converts both structures into
cleaned single protein chains, optionally adds hydrogens using OpenMM, remaps
the input amino acid substitution(s) to the cleaned residue number(s),
and outputs cleaned structure files together with a JSON report summarizing
residue mapping and cleaning statistics for both wild-type and mutant
structures.

'''

# example usage:
'''
Example command:

enzywizard-mut-clean -w examples/input/1ZG4_WT.cif -m examples/input/1ZG6_S70G.cif -s S70G -o examples/output/

'''

# input parameters:
'''
-w, --wt_input_path
Required.
Path to the input wild-type protein structure file in CIF or PDB format.

-m, --mut_input_path
Required.
Path to the input mutant protein structure file in CIF or PDB format.

-s, --mutation
Required.
Amino acid substitution(s) describing the mutation relationship between the
wild-type and mutant structures.
Examples:
A123V
A123V,G456D
Multiple mutations should be separated by ','.

-o, --output_dir
Required.
Path to the output directory for saving cleaned structure files and the JSON report.

--no_add_H
Optional flag.
Disable hydrogen addition using OpenMM.
By default, hydrogens are added.

--pH
Optional.
pH value used for hydrogen addition.
Default: 7.0
Valid range: 0.0 to 14.0
'''

# output content:
'''
The program outputs the following files into the output directory:

1. Cleaned structure files for the wild-type structure
   - cleaned_{wt_name}.cif
   - cleaned_{wt_name}.pdb
   - cleaned_{wt_name}.fasta

2. Cleaned structure files for the mutant structure
   - cleaned_{mut_name}.cif
   - cleaned_{mut_name}.pdb
   - cleaned_{mut_name}.fasta

3. A JSON report
   - mut_clean_report_{wt_name}_to_{mut_name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_mut_clean"

   - "amino_acid_substitution"
     The original amino acid substitution(s) provided by the user.

   - "cleaned_amino_acid_substitution"
     The amino acid substitution(s) after remapping from the original residue
     numbering system to the cleaned residue numbering system.

   - "wt_amino_acid_mapping_old_to_new"
     A list describing how residues in the original wild-type structure
     correspond to residues in the cleaned wild-type structure.

     Each entry contains:
     - "old_residue"
       Information for the original residue before cleaning:
       - "aa_id": original residue index
       - "aa_name": original residue one-letter amino acid code
       - "hydrogen_atom_count": number of hydrogen atoms found in the original residue

     - "new_residue"
       Information for the cleaned residue after cleaning:
       - "aa_id": new residue index after continuous renumbering
       - "aa_name": cleaned residue one-letter amino acid code
       - "hydrogen_atom_count": number of hydrogen atoms in the cleaned residue

   - "wt_clean_statistics"
     A dictionary summarizing the overall cleaning process for the wild-type structure.

     It includes:
     - "changed_resname"
       Number of residues whose names were standardized using the MODRES mapping.

     - "removed_nonstd"
       Number of non-standard residues removed.

     - "removed_missing_bb"
       Number of residues removed because required backbone atoms
       (N, CA, C) were missing.

     - "removed_missing_heavy_atoms"
       Number of residues removed because one or more required heavy atoms
       were missing.

     - "removed_unexpected_heavy_atoms"
       Number of residues removed because unexpected heavy atoms were present.

     - "removed_bad_occ"
       Number of residues removed because selected atoms had invalid occupancy.

     - "removed_inscodes"
       Number of residues whose original insertion codes were present and then removed
       during residue renumbering.

     - "kept_residues"
       Number of residues kept in the final cleaned wild-type structure.

   - "mut_amino_acid_mapping_old_to_new"
     A list describing how residues in the original mutant structure
     correspond to residues in the cleaned mutant structure.

     Each entry contains:
     - "old_residue"
       Information for the original residue before cleaning:
       - "aa_id": original residue index
       - "aa_name": original residue one-letter amino acid code
       - "hydrogen_atom_count": number of hydrogen atoms found in the original residue

     - "new_residue"
       Information for the cleaned residue after cleaning:
       - "aa_id": new residue index after continuous renumbering
       - "aa_name": cleaned residue one-letter amino acid code
       - "hydrogen_atom_count": number of hydrogen atoms in the cleaned residue

   - "mut_clean_statistics"
     A dictionary summarizing the overall cleaning process for the mutant structure.

     It includes:
     - "changed_resname"
       Number of residues whose names were standardized using the MODRES mapping.

     - "removed_nonstd"
       Number of non-standard residues removed.

     - "removed_missing_bb"
       Number of residues removed because required backbone atoms
       (N, CA, C) were missing.

     - "removed_missing_heavy_atoms"
       Number of residues removed because one or more required heavy atoms
       were missing.

     - "removed_unexpected_heavy_atoms"
       Number of residues removed because unexpected heavy atoms were present.

     - "removed_bad_occ"
       Number of residues removed because selected atoms had invalid occupancy.

     - "removed_inscodes"
       Number of residues whose original insertion codes were present and then removed
       during residue renumbering.

     - "kept_residues"
       Number of residues kept in the final cleaned mutant structure.

   This report helps users track:
   - how residues in both structures were retained and renumbered,
   - how residue names were standardized,
   - how hydrogen content changed after optional hydrogen addition,
   - and how the input mutation positions were translated into the cleaned
     residue numbering system.
'''

# Process:
'''
This command processes the input wild-type and mutant protein structures as follows:

1. Load the input structures
   - Read the wild-type CIF or PDB file using Biopython (Bio.PDB).
   - Read the mutant CIF or PDB file using Biopython (Bio.PDB).
   - Resolve the protein names from the input filenames.

2. Validate basic input conditions
   - Check that both input files exist.
   - Check that the wild-type and mutant filenames are valid and not identical.
   - Extract a single chain from each structure.
   - Check that the wild-type and mutant sequence lengths are valid.
   - Check that the input amino acid substitution format is valid.
   - Check that mutation positions fall within the residue index ranges of both structures.
   - Check that the pH value is within the valid range.

3. Clean both structures independently (Biopython-based processing)
   - Extract a single chain from each structure using Biopython structure utilities.
   - Standardize residue names using the MODRES mapping.
   - Remove residues with missing backbone atoms (N, CA, C).
   - Remove residues with missing required heavy atoms.
   - Remove residues with unexpected heavy atoms.
   - Remove residues with invalid occupancy.
   - Remove insertion codes.
   - Repair discontinuous residue numbering by rebuilding residue indices.
   - Renumber all kept residues continuously starting from 1.
   - Rebuild each output structure as a single chain with chain ID A.

4. Remap amino acid substitution
   - Build residue mapping relationships from original residue numbering to cleaned
     residue numbering for both wild-type and mutant structures.
   - Convert the input amino acid substitution(s) into cleaned amino acid
     substitution(s) using the residue mapping results.

5. Optionally add hydrogens
   - Convert each cleaned Biopython structure into an OpenMM object.
   - Use OpenMM Modeller.addHydrogens() with a specified pH.
   - Convert the hydrogen-added structures back into Biopython structures.

6. Validate the cleaned structures
   - Check that each cleaned structure contains exactly one model and one chain.
   - Check that the cleaned chain ID is A.
   - Check that residue numbering is continuous.
   - Check that only standardized standard amino acid residues remain.
   - Check that required backbone and heavy atoms are present.
   - Check that no unexpected heavy atoms remain.
   - Verify that mapped residue coordinates are preserved after cleaning.
   - Check that the cleaned wild-type and mutant sequence lengths are equal.

7. Save outputs
   - Save the cleaned wild-type structure in CIF format (Biopython MMCIFIO).
   - Save the cleaned wild-type structure in PDB format (Biopython PDBIO).
   - Extract and save the cleaned wild-type amino acid sequence in FASTA format.
   - Save the cleaned mutant structure in CIF format (Biopython MMCIFIO).
   - Save the cleaned mutant structure in PDB format (Biopython PDBIO).
   - Extract and save the cleaned mutant amino acid sequence in FASTA format.
   - Generate and save the JSON report summarizing mutation remapping,
     residue mappings, and cleaning statistics for both structures.
'''

# dependencies:
'''
- Biopython
- OpenMM
- NumPy
'''

# references:
'''
- Biopython:
  https://biopython.org/

- OpenMM:
  https://openmm.org/

- wwPDB Chemical Component Dictionary / MODRES-related residue standardization resource:
  https://www.wwpdb.org/data/ccd

- Rosetta structure preparation overview:
  https://docs.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-structures
'''


