from __future__ import annotations
from argparse import Namespace
from ..services.hydrocluster_service import run_hydrocluster_service


def add_hydrocluster_parser(subparsers) -> None:
    parser = subparsers.add_parser("hydrocluster",help="EnzyWizard-HydroCluster: Identify hydrophobic clusters from input CIF/PDB file and generate a detailed JSON report.")
    parser.add_argument("-i","--input_path", required=True, help="Path to the input cleaned protein structure file in CIF or PDB format.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to the output directory for saving the JSON report.")
    parser.add_argument("--cutoff", type=float, default=10.0, help="Minimum residue-residue contact area cutoff for hydrophobic cluster connection.")
    parser.set_defaults(func=run_hydrocluster)

def run_hydrocluster(args: Namespace) -> None:
    run_hydrocluster_service(input_path=args.input_path, output_dir=args.output_dir, cutoff_area=args.cutoff)

# ==============================
# Command: enzywizard-hydrocluster
# ==============================

# brief introduction:
'''
EnzyWizard-HydroCluster is a command-line tool for identifying hydrophobic
clusters from a cleaned protein structure and generating a detailed JSON report.
It detects residue-residue hydrophobic contacts among ILE, VAL, and LEU residues
based on estimated side-chain surface-contact areas, groups connected residues
into hydrophobic clusters, and computes overall statistics summarizing the number
and area of detected clusters across the protein.

'''

# example usage:
'''
Example command:

enzywizard-hydrocluster -i examples/input/cleaned_3GP6.cif -o examples/output/

'''

# input parameters:
'''
-i, --input_path
Required.
Path to the input cleaned protein structure file in CIF or PDB format.

-o, --output_dir
Required.
Path to the output directory for saving the JSON report.

--cutoff
Optional.
Minimum residue-residue contact area cutoff for hydrophobic cluster connection.
Default: 10.0
Must be a positive number.
'''

# output content:
'''
The program outputs the following file into the output directory:

1. A JSON report
   - hydrocluster_report_{name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_hydrocluster"

   - "hydrophobic_cluster_statistics"
     A dictionary summarizing hydrophobic cluster statistics over the full protein.

     It includes:
     - "cluster_num"
       Total number of detected hydrophobic clusters.

     - "max_cluster_area"
       Maximum total cluster area among all detected hydrophobic clusters.

     - "total_cluster_area"
       Sum of total areas over all detected hydrophobic clusters.

   - "hydrophobic_cluster"
     A list describing hydrophobic clusters detected in the cleaned protein
     structure.

     Each entry contains:
     - "area"
       Total estimated contact area of the hydrophobic cluster.

     - "residues"
       A list of residues involved in the hydrophobic cluster.

       Each residue entry contains:
       - "aa_id"
         Residue index in the cleaned structure.

       - "aa_name"
         Residue one-letter amino acid code.
'''

# Process:
'''
This command processes the input cleaned protein structure as follows:

1. Load the input structure
   - Read the cleaned CIF or PDB file using Biopython (Bio.PDB).
   - Resolve the protein name from the input filename.

2. Validate basic input conditions
   - Check that the input file exists.
   - Check that the cutoff value is a positive number.
   - Validate that the input structure satisfies the cleaned-structure requirement.

3. Extract residues and atoms for hydrophobic cluster calculation
   - Extract the single chain from the cleaned structure.
   - Identify all ILE, VAL, and LEU residues in the chain.
   - Extract side-chain non-hydrogen atoms from these ILE/VAL/LEU residues.
   - Extract all protein non-hydrogen atoms in the chain as possible neighboring atoms.

4. Estimate residue-residue hydrophobic contact areas
   - For each ILE/VAL/LEU side-chain atom, generate evenly distributed sample
     points on a sphere centered at the atom.
   - Use Biopython NeighborSearch to identify nearby non-hydrogen atoms.
   - Determine which sphere surface points are covered by neighboring atoms.
   - Assign each covered point to one neighboring atom and estimate atom-level
     contact areas from the number of covered sample points.
   - Keep only contacts that map to ILE/VAL/LEU side-chain atoms and accumulate
     them into a residue-level contact area matrix.

5. Build hydrophobic clusters
   - Build a residue contact graph in which each node represents one ILE, VAL,
     or LEU residue.
   - Add directed residue-residue edges when the estimated contact area is above
     the defined cutoff.
   - Identify hydrophobic clusters as weakly connected components in the residue
     contact graph.
   - Compute the total area of each hydrophobic cluster and sort clusters by area.

6. Compute summary statistics
   - Count the total number of detected hydrophobic clusters.
   - Calculate the maximum hydrophobic cluster area.
   - Calculate the total hydrophobic cluster area across the protein.

7. Save outputs
   - Generate and save a JSON report containing both hydrophobic cluster details
     and overall hydrophobic cluster statistics.
'''

# dependencies:
'''
- Biopython
- NumPy
- SciPy
- NetworkX
'''

# references:
'''
- Biopython:
  https://biopython.org/

- SciPy:
  https://scipy.org/

- NetworkX:
  https://networkx.org/

- Protlego hydrophobic cluster implementation:
  https://github.com/Hoecker-Lab/protlego/blob/master/protlego/structural/clusters.py

- Protlego:
  https://github.com/Hoecker-Lab/protlego
'''