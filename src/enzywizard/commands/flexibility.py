from __future__ import annotations
from argparse import Namespace
from ..services.flexibility_service import run_flexibility_service


def add_flexibility_parser(subparsers) -> None:
    parser = subparsers.add_parser("flexibility", help="Calculate protein flexibility (RMSF) from input CIF/PDB file.")
    parser.add_argument("-i", "--input_path", required=True, help="Path to input CIF/PDB file.")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to a directory for outputting a JSON report.")
    parser.add_argument("--method",type=str,choices=["ANM", "GNM"],default="ANM",help="Method for RMSF calculation: ANM or GNM (default: ANM).")
    parser.add_argument("--cutoff",type=float,default=15.0,help="Distance cutoff used to determine the residue connection in ProDy (default: 15.0).")
    parser.add_argument("--n_modes",type=int,default=20,help="Number of low-frequency normal modes used for RMSF calculation (default: 20).")

    parser.set_defaults(func=run_flexibility)

def run_flexibility(args: Namespace) -> None:
    run_flexibility_service(input_path=args.input_path,output_dir=args.output_dir,cutoff=args.cutoff,n_modes=args.n_modes,method=args.method)

# input parameters:
'''
-i --input_path Required. Input cleaned CIF/PDB protein structure file.

-o --output_dir Required. Output directory for saving a JSON report.

--method Optional. Method for RMSF calculation (default: ANM).
Supported values:
- ANM: Anisotropic Network Model
- GNM: Gaussian Network Model

--cutoff Optional. Distance cutoff for building the residue connection in ProDy (default: 15.0).
Residues whose CA atoms are within this cutoff are considered connected.
This parameter controls the connectivity density of the elastic network.
A smaller cutoff gives a sparser network, while a larger cutoff gives a denser network.

--n_modes Optional. Number of low-frequency normal modes used for RMSF calculation (default: 20).
These modes represent collective motions of the protein.
Using more modes includes more motion information, while using fewer modes
focuses more on the largest-scale global motions.
'''

# output content:
'''
The program outputs a JSON report that records:

1. "output_type": "enzywizard_flexibility"

2. "protein_rmsf": a list of per-residue RMSF records, each including:
   - aa_id: residue sequence index
   - aa_name: residue name
   - rmsf: calculated RMSF value
'''

# functionality/process
'''
It processes a cleaned protein structure by:

1. Loading the input structure in Biopython format;

2. Checking whether the input structure satisfies the EnzyWizard cleaned
   structure requirements, including:
   - exactly one model and one chain with chain ID "A"
   - continuous residue numbering starting from 1
   - no insertion codes
   - only standard amino acid residues
   - complete backbone atoms (N, CA, C)
   - complete required heavy atoms
   - valid occupancy values

3. Extracting CA coordinates for all amino acid residues;

4. Building an elastic network model in ProDy based on the selected method;

    ANM (Anisotropic Network Model): 

     Builds a Hessian matrix from the CA-based elastic network, solves the
     low-frequency normal modes, computes residue square fluctuations from
     these modes, and converts them to RMSF.

     ANM models residue motions as 3D directional fluctuations in an elastic
     network. It preserves motion direction information, so it is more suitable
     for describing anisotropic and collective structural motions.
    
    GNM (Gaussian Network Model)

     Builds a Kirchhoff matrix from the CA-based elastic network, solves the
     low-frequency normal modes, computes residue square fluctuations from
     these modes, and converts them to RMSF.

     GNM models residue motions as isotropic fluctuations in a residue contact
     network. It does not describe motion directions, but captures the relative
     mobility of residues efficiently and robustly.
    

5. Calculating square fluctuations from the selected normal modes and converting
   them to RMSF values;

6. Generating a structured JSON report containing per-residue RMSF values.
'''

# dependency:
'''
Biopython
ProDy
NumPy
'''

# reference:
'''
ProDy documentation
https://prody.csb.pitt.edu/
'''

