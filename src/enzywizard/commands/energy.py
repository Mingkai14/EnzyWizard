from __future__ import annotations
from argparse import Namespace
from ..services.energy_service import run_energy_service


def add_energy_parser(subparsers) -> None:
    parser = subparsers.add_parser("energy",help="Calculate energy terms from input CIF/PDB file.")
    parser.add_argument("-i","--input_path", required=True, help="Path to input CIF/PDB file.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to a directory for outputting a JSON report.")
    parser.add_argument("--minimize_energy",type=lambda x: str(x).lower() in ["true", "1", "yes"],default=True,help="Whether to perform an energy minimization before energy evaluation (True/False, default: True).")
    parser.add_argument("--minimization_iteration",type=int,default=1000,help="Maximum number of iterations for energy minimization (default: 1000). A smaller value may result in energy values closer to the instantaneous (non-minimized) state.")
    parser.set_defaults(func=run_energy)


def run_energy(args: Namespace) -> None:
    run_energy_service(input_path=args.input_path,output_dir=args.output_dir,minimize_energy=args.minimize_energy,minimization_iteration=args.minimization_iteration)

# input parameters:
'''
-i --input_path Required. Input cleaned CIF/PDB protein structure file.

-o --output_dir Required. Output directory for saving a JSON report.

--minimize_energy Optional. Whether to perform an energy minimization before energy evaluation (default: True).

--minimization_iteration Optional. Maximum number of iterations for energy minimization (default: 1000). A smaller value may result in energy values closer to the instantaneous (non-minimized) state.
'''

# output content:
'''
The program outputs a JSON report that records:

1. "output_type": "enzywizard_energy"

2. "energy_terms": a dictionary of energy terms calculated by OpenMM,
   including:
   - total_potential_energy
   - harmonic_bond_force
   - harmonic_angle_force
   - custom_bond_force
   - custom_torsion_force
   - custom_nonbonded_force
   - nonbonded_force
   - periodic_torsion_force
   - cmap_torsion_force
'''
# functionality/process
'''
It processes a cleaned protein structure by:

1. Loading the input structure in both Biopython and OpenMM formats;

2. Checking whether the input structure satisfies the EnzyWizard cleaned
   structure requirements, including:
   - exactly one model and one chain with chain ID "A"
   - continuous residue numbering starting from 1
   - no insertion codes
   - only standard amino acid residues
   - complete backbone atoms (N, CA, C)
   - complete required heavy atoms
   - valid occupancy values

3. Building an OpenMM system with the specified force field;

4. Optionally performing an energy minimization before energy evaluation;

5. Calculating the total potential energy and individual force-field energy
   terms from the OpenMM context;

6. Generating a structured JSON report containing the calculated energy terms.
'''

# dependency:
'''
Biopython
OpenMM
'''

# reference:
'''
OpenMM documentation
https://openmm.org/
'''