from __future__ import annotations
from argparse import Namespace
from ..services.energy_service import run_energy_service


def add_energy_parser(subparsers) -> None:
    parser = subparsers.add_parser("energy",help="Calculate molecular mechanics energy terms from a cleaned protein structure and generate a detailed JSON report")
    parser.add_argument("-i","--input_path", required=True, help="Path to the input cleaned protein structure file in CIF or PDB format.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to the output directory for saving the JSON report.")
    parser.add_argument("--not_minimize_energy",action="store_false",dest="minimize_energy",help="Disable performing an energy minimization before energy evaluation (default: enabled).")
    parser.set_defaults(minimize_energy=True)
    parser.add_argument("--minimization_iteration",type=int,default=1000,help="Maximum number of iterations for energy minimization (default: 1000). A smaller value may result in energy values closer to the instantaneous (non-minimized) state.")
    parser.set_defaults(func=run_energy)


def run_energy(args: Namespace) -> None:
    run_energy_service(input_path=args.input_path,output_dir=args.output_dir,minimize_energy=args.minimize_energy,minimization_iteration=args.minimization_iteration)


# ==============================
# Command: enzywizard-energy
# ==============================

# brief introduction:
'''
EnzyWizard-Energy is a command-line tool for calculating molecular mechanics
energy terms from a cleaned protein structure and generating a detailed JSON report.
It evaluates the total potential energy and multiple force-field energy components
using OpenMM. Optionally, it can perform an energy minimization before energy
evaluation, so the reported energy values can reflect either the minimized structure
state or a state closer to the original input conformation.

'''

# example usage:
'''
Example command:

enzywizard-energy -i examples/input/cleaned_3GP6.cif -o examples/output/

'''

# input parameters:
'''
-i, --input_path
Required.
Path to the input cleaned protein structure file in CIF or PDB format.

-o, --output_dir
Required.
Path to the output directory for saving the JSON report.

--not_minimize_energy
Optional.
Disable performing an energy minimization before energy evaluation.
By default, energy minimization is enabled.

--minimization_iteration
Optional.
Maximum number of iterations for energy minimization.
Default: 1000.
A smaller value may result in energy values closer to the instantaneous
(non-minimized) structure state.
'''

# output content:
'''
The program outputs the following file into the output directory:

1. A JSON report
   - energy_report_{name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_energy"

   - "energy_terms"
     A dictionary containing molecular mechanics energy terms calculated
     from the cleaned protein structure using OpenMM.

     It includes:
     - "total_potential_energy"
       Total potential energy of the structure.

     - "harmonic_bond_force"
       Energy contribution from harmonic bond terms.

     - "harmonic_angle_force"
       Energy contribution from harmonic angle terms.

     - "custom_bond_force"
       Energy contribution from OpenMM custom bond terms, if present in the system.

     - "custom_torsion_force"
       Energy contribution from OpenMM custom torsion terms, if present in the system.

     - "custom_nonbonded_force"
       Energy contribution from OpenMM custom nonbonded interaction terms, if present
       in the system.

     - "nonbonded_force"
       Energy contribution from standard nonbonded interactions.

     - "periodic_torsion_force"
       Energy contribution from periodic torsion terms.

     - "cmap_torsion_force"
       Energy contribution from CMAP torsion correction terms, if present
       in the system.

   All energy values are reported in kilojoule per mole (kJ/mol).
'''

# Process:
'''
This command processes the input cleaned protein structure as follows:

1. Load the input structure
   - Read the cleaned CIF or PDB file using Biopython (Bio.PDB).
   - Load the same structure into OpenMM format for force-field-based
     energy evaluation.
   - Resolve the protein name from the input filename.

2. Validate basic input conditions
   - Check that the input file exists.
   - Validate that the input structure satisfies the cleaned-structure requirement.

3. Build the OpenMM system
   - Load the specified OpenMM force field.
   - Construct an OpenMM Modeller object from the cleaned structure.
   - Create an OpenMM system using the protein topology and coordinates.
   - Assign force groups to individual force terms.

4. Initialize the OpenMM context
   - Create the OpenMM integrator and CPU platform context.
   - Set the cleaned protein coordinates into the context.

5. Optionally perform energy minimization
   - If minimization is enabled, run OpenMM local energy minimization
     before evaluating the energy terms.
   - The maximum number of minimization iterations is controlled by
     the user-defined parameter.

6. Calculate energy terms
   - Evaluate the total potential energy of the structure.
   - Evaluate individual force-field energy terms by querying each
     assigned OpenMM force group separately.

7. Save outputs
   - Generate and save a JSON report containing the calculated
     energy terms.
'''

# dependencies:
'''
- Biopython
- OpenMM
'''

# references:
'''
- OpenMM:
  https://openmm.org/

- OpenMM documentation:
  https://docs.openmm.org/
'''