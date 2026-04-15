from __future__ import annotations
from argparse import Namespace
from ..services.clean_service import run_clean_service

def add_clean_parser(subparsers) -> None:
    parser = subparsers.add_parser("clean",help="Clean a CIF/PDB structure file.")
    parser.add_argument("-i","--input_path", required=True, help="Path to input CIF/PDB file.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to a directory for outputting cleaned CIF, PDB, and FASTA files and a JSON report.")
    parser.add_argument("--no_add_H",action="store_false",dest="add_H",help="Disable adding hydrogens using OpenMM (default: enabled).")
    parser.set_defaults(add_H=True)
    parser.add_argument("--pH",type=float,default=7.0,help="pH value for hydrogen addition (default: 7.0).")
    parser.set_defaults(func=run_clean)

def run_clean(args: Namespace) -> None:
    run_clean_service(input_path=args.input_path, output_dir=args.output_dir, add_H=args.add_H, pH=args.pH)


# ==============================
# Command: enzywizard-clean
# ==============================

# brief introduction:
'''
EnzyWizard-Clean is a command-line tool for cleaning an input protein structure
file in CIF or PDB format. It standardizes residue names, removes problematic
residues (non-standard residues, residues with missing backbone atoms,
residues with incomplete heavy atoms, residues with invalid occupancy), 
repairs residue order by renumbering residues continuously, converts the 
structure into a cleaned protein chain, optionally
adds hydrogens using OpenMM, and outputs cleaned structure files together with
a JSON report summarizing residue mapping and cleaning statistics.

'''

# example usage:
'''
Example command:

enzywizard-clean -i examples/input/3GP6.cif -o examples/output/

'''

# input parameters:
'''
-i, --input_path
Required.
Path to the input protein structure file in CIF or PDB format.

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

1. Cleaned structure files
   - cleaned_{name}.cif
   - cleaned_{name}.pdb
   - cleaned_{name}.fasta

2. A JSON report
   - clean_report_{name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_clean"

   - "amino_acid_mapping_old_to_new"
     A list describing how residues in the original structure correspond to
     residues in the cleaned structure.

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

     This mapping helps users track:
     - which residues were kept,
     - how residue numbering changed,
     - whether residue names were standardized,
     - and how hydrogen content changed after optional hydrogen addition.

   - "clean_statistics"
     A dictionary summarizing the overall cleaning process.

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

     - "removed_bad_occ"
       Number of residues removed because selected key atoms had invalid occupancy.

     - "removed_inscodes"
       Number of residues whose original insertion codes were present and then removed
       during residue renumbering.

     - "kept_residues"
       Number of residues kept in the final cleaned structure.
'''

# Process:
'''
This command processes the input protein structure as follows:

1. Load the input structure
   - Read the CIF or PDB file using Biopython (Bio.PDB).
   - Resolve the protein name from the input filename.

2. Validate basic input conditions
   - Check that the input file exists.
   - Check that the pH value is within the valid range.

3. Clean the structure (Biopython-based processing)
   - Extract a single chain using Biopython structure utilities.
   - Standardize non-standard residues using the internal MODRES mapping (wwPDB Chemical Component Dictionary).
   - Remove residues with missing backbone atoms (N, CA, C).
   - Remove residues with missing required heavy atoms.
   - Remove residues with invalid occupancy values.
   - Remove insertion codes.
   - Repair discontinuous residue numbering by rebuilding residue indices.
   - Renumber all kept residues continuously starting from 1.

4. Optionally add hydrogens
   - Convert the cleaned Biopython structure into an OpenMM object.
   - Use OpenMM Modeller.addHydrogens() with a specified pH.
   - Convert the hydrogen-added structure back into a Biopython structure.

5. Validate the cleaned structure

6. Save outputs
   - Save the cleaned structure in CIF format (Biopython MMCIFIO).
   - Save the cleaned structure in PDB format (Biopython PDBIO).
   - Extract and save the cleaned amino acid sequence in FASTA format.
   - Generate and save the JSON report summarizing residue mapping and statistics.
'''

# dependencies:
'''
- Biopython
- OpenMM
- NumPy
- MODRES residue-name mapping table
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