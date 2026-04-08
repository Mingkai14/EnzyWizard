from __future__ import annotations
from argparse import Namespace
from ..services.clean_service import run_clean_service

def add_clean_parser(subparsers) -> None:
    parser = subparsers.add_parser("clean",help="Clean a CIF/PDB structure file.")
    parser.add_argument("-i","--input_path", required=True, help="Path to input CIF/PBD file.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to a directory for outputting cleaned CIF, PDB, and FASTA files and a JSON report.")
    parser.add_argument( "--add_H", type=lambda x: str(x).lower() in ["true", "1", "yes"], default=True, help="Whether to add hydrogens using OpenMM (True/False, default: True)." )
    parser.add_argument("--pH",type=float,default=7.0,help="pH value for hydrogen addition (default: 7.0).")

    parser.set_defaults(func=run_clean)

def run_clean(args: Namespace) -> None:
    run_clean_service(input_path=args.input_path, output_dir=args.output_dir, add_H=args.add_H, pH=args.pH)


# =========================
# Command: enzywizard clean
# =========================

# input parameters:
'''
-i --input_path Required. Path to input protein structure file (CIF or PDB).

-o --output_dir Required. Directory to save cleaned structure and sequence files and report.

--add_H Optional. Whether to add hydrogens using OpenMM (True/False). Default: True.

--pH Optional. pH value used for hydrogen addition. Default: 7.0.
'''

# output content:
'''
The program outputs:

1. Cleaned structure files:
   - cleaned_{name}.cif
   - cleaned_{name}.pdb
   - cleaned_{name}.fasta

2. A JSON report:
   - clean_report_{name}.json

The JSON report contains:

- "output_type": "enzywizard_clean"

- "amino_acid_mapping_old_to_new":
  A list of mappings from original residues to cleaned residues.
  Each entry includes:
    - old_residue:
        - aa_id
        - aa_name
        - hydrogen_atom_count
    - new_residue:
        - aa_id
        - aa_name
        - hydrogen_atom_count

- "clean_statistics":
  Statistics of the cleaning process, including:
    - changed_resname
    - removed_nonstd
    - removed_missing_bb
    - removed_missing_heavy_atoms
    - removed_bad_occ
    - removed_inscodes
    - kept_residues
'''

# functionality / process:
'''
This command processes a protein structure as follows:

1. Load structure using Biopython.

2. Clean the structure:
   - Extract a single chain.
   - Standardize residue names using MODRES mapping.
   - Remove non-standard residues.
   - Remove residues with missing backbone atoms (N, CA, C).
   - Remove residues with incomplete heavy atoms.
   - Remove residues with invalid occupancy.
   - Remove insertion codes.
   - Renumber residues consecutively.

3. (Optional) Add hydrogens:
   - Convert structure to OpenMM PDBFile.
   - Add hydrogens using Modeller.addHydrogens() with specified pH and force field.
   - Convert back to Biopython Structure.

4. Validate structure:
   - Check structural integrity.
   - Verify that residue coordinates (CA atoms) are unchanged after cleaning.

5. Save outputs:
   - Cleaned CIF, PDB, and FASTA files.
   - JSON report with mapping and statistics.
'''

# dependencies:
'''
- Biopython
- OpenMM
- MODRES mapping table (https://www.wwpdb.org/data/ccd)
'''

# reference:
'''
- https://docs.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-structures
'''