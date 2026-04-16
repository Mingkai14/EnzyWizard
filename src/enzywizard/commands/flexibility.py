from __future__ import annotations
from argparse import Namespace
from ..services.flexibility_service import run_flexibility_service


def add_flexibility_parser(subparsers) -> None:
    parser = subparsers.add_parser("flexibility", help="Estimate residue flexibility from a cleaned protein structure and generate a detailed JSON report.")
    parser.add_argument("-i", "--input_path", required=True, help="Path to the input cleaned protein structure file in CIF or PDB format.")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output directory for saving the JSON report.")
    parser.add_argument("--method",type=str,choices=["ANM", "GNM"],default="ANM",help="Method for RMSF calculation: ANM or GNM (default: ANM).")
    parser.add_argument("--cutoff",type=float,default=15.0,help="Distance cutoff used to determine the residue connection in ProDy (default: 15.0).")
    parser.add_argument("--n_modes",type=int,default=20,help="Number of low-frequency normal modes used for RMSF calculation (default: 20).")

    parser.set_defaults(func=run_flexibility)

def run_flexibility(args: Namespace) -> None:
    run_flexibility_service(input_path=args.input_path,output_dir=args.output_dir,cutoff=args.cutoff,n_modes=args.n_modes,method=args.method)

# ==============================
# Command: enzywizard-flexibility
# ==============================

# brief introduction:
'''
EnzyWizard-Flexibility is a command-line tool for estimating residue
flexibility from a cleaned protein structure and generating a detailed JSON report.
It computes per-residue root mean square fluctuation (RMSF) values using elastic
network models implemented in ProDy, including the Anisotropic Network Model (ANM)
and Gaussian Network Model (GNM). These models capture intrinsic protein dynamics
based on the topology of the structure and provide insight into residue mobility
and collective motions.

'''

# example usage:
'''
Example command:

enzywizard-flexibility -i examples/input/cleaned_3GP6.cif -o examples/output/

'''

# input parameters:
'''
-i, --input_path
Required.
Path to the input cleaned protein structure file in CIF or PDB format.

-o, --output_dir
Required.
Path to the output directory for saving the JSON report.

--method
Optional.
Method for RMSF calculation (default: ANM).
Supported values:
- ANM: Anisotropic Network Model
- GNM: Gaussian Network Model

--cutoff
Optional.
Distance cutoff for building the residue connection in ProDy (default: 15.0).
Residues whose CA atoms are within this cutoff are considered connected.
This parameter controls the connectivity density of the elastic network.
A smaller cutoff gives a sparser network, while a larger cutoff gives a denser network.

--n_modes
Optional.
Number of low-frequency normal modes used for RMSF calculation (default: 20).
These modes represent collective motions of the protein.
Using more modes includes more motion information, while using fewer modes
focuses more on the largest-scale global motions.
'''

# output content:
'''
The program outputs the following file into the output directory:

1. A JSON report
   - flexibility_report_{name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_flexibility"

   - "protein_rmsf"
     A list describing residue-level flexibility for each residue in the
     cleaned protein structure.

     Each entry contains:
     - "aa_id"
       Residue index in the cleaned structure.

     - "aa_name"
       Residue one-letter amino acid code.

     - "rmsf"
       Root mean square fluctuation (RMSF) value estimated from the selected
       elastic network model.
'''

# Process:
'''
This command processes the input cleaned protein structure as follows:

1. Load the input structure
   - Read the cleaned CIF or PDB file using Biopython (Bio.PDB).
   - Resolve the protein name from the input filename.

2. Validate basic input conditions
   - Check that the input file exists.
   - Validate that the input structure satisfies the cleaned-structure requirement.

3. Extract structural information
   - Extract the single chain from the cleaned structure.
   - Retrieve all residues in chain order.
   - Extract CA atom coordinates for each residue.
   - Ensure sufficient residues are available for flexibility calculation.

4. Build elastic network model
   - Construct an elastic network model using ProDy based on the selected method:

     ANM (Anisotropic Network Model):
     - Builds a Hessian matrix from the CA-based elastic network.
     - Solves low-frequency normal modes.
     - Computes square fluctuations from these modes.
     - Converts square fluctuations into RMSF values.
     - Captures directional (anisotropic) motions of residues.

     GNM (Gaussian Network Model):
     - Builds a Kirchhoff matrix from the CA-based elastic network.
     - Solves low-frequency normal modes.
     - Computes square fluctuations from these modes.
     - Converts square fluctuations into RMSF values.
     - Captures isotropic residue mobility without directional information.

5. Calculate RMSF values
   - Convert square fluctuations into RMSF values by taking square roots.
   - Ensure consistency between residue count and RMSF results.

6. Assemble residue-level results
   - For each residue, record:
     - residue index (aa_id)
     - amino acid type (aa_name)
     - RMSF value (rmsf)

7. Save outputs
   - Generate and save a JSON report containing per-residue flexibility values.
'''

# dependencies:
'''
- Biopython
- ProDy
- NumPy
'''

# references:
'''
- ProDy:
  https://prody.csb.pitt.edu/

- ProDy documentation:
  https://prody.csb.pitt.edu/manual/

- Elastic Network Models:
  Bahar et al., "Direct evaluation of thermal fluctuations in proteins using a single-parameter harmonic potential", Folding & Design (1997)
'''

