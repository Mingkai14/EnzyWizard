from __future__ import annotations
from argparse import Namespace
from ..services.hydrocluster_service import run_hydrocluster_service


def add_hydrocluster_parser(subparsers) -> None:
    parser = subparsers.add_parser("hydrocluster",help="Calculate hydrophobic clusters from input CIF/PDB file.")
    parser.add_argument("-i","--input_path", required=True, help="Path to input CIF/PDB file.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to a directory for outputting a JSON report.")
    parser.add_argument("-c", "--cutoff", type=float, default=10.0, help="Minimum contact area cutoff for hydrophobic cluster residue-residue connection (default: 10.0).")
    parser.set_defaults(func=run_hydrocluster)

def run_hydrocluster(args: Namespace) -> None:
    run_hydrocluster_service(input_path=args.input_path, output_dir=args.output_dir, cutoff_area=args.cutoff)

# input parameters:
'''
-i --input_path required input cleaned CIF/PDB protein structure file;
-o --output_dir required output directory to save JSON report
-c --cutoff optional minimum contact area cutoff for hydrophobic cluster residue-residue connection (default: 10.0)
'''

# output content:
'''
The program outputs a JSON report that records:

1. "output_type": enzywizard_hydrocluster,

2. "hydrophobic_cluster": a list of hydrophobic clusters, where each cluster includes:
   - cluster surface-contact area (area)
   - a list of residues involved in the cluster (residues), and for each residue:
     - amino acid index (aa_id)
     - amino acid name (aa_name)
'''

# functionality/process
'''
It processes a protein structure by:

1. Extracting a single protein chain from the input structure;

2. Identifying all hydrophobic residues of type ILE, VAL, and LEU in the chain,
   and extracting their side-chain non-hydrogen atoms, which are used as the
   basic units for hydrophobic cluster calculation;

3. Extracting all non-hydrogen atoms in the protein chain as possible neighbors
   around each ILE/VAL/LEU side-chain atom;

4. For each ILE/VAL/LEU side-chain atom:
   - treating the atom as a sphere with carbon atomic radius plus solvent probe
     radius,
   - generating evenly distributed sample points on the sphere surface,
   - searching all nearby non-hydrogen atoms within a fixed cutoff distance,
   - checking, for each sampled surface point, whether it falls inside the
     sphere of any neighboring atom;

5. If one sampled surface point is covered by multiple neighboring atoms,
   assigning that point to the nearest neighboring atom center, so that each
   covered point contributes to only one atom-level contact;

6. Estimating atom-level contact area by counting how many sampled surface points
   are covered by each neighboring atom, and multiplying the number of covered
   points by the surface area represented by one sample point;

7. Keeping only atom-level contacts in which both the source atom and the
   contacting neighboring atom belong to side-chain non-hydrogen atoms of
   ILE/VAL/LEU residues, and summing these atom-level contact areas into a
   residue-level contact area matrix;

8. Building a directed residue contact graph, in which:
   - each node represents one ILE/VAL/LEU residue,
   - each directed edge represents a residue-to-residue contact whose accumulated
     contact area is above a defined cutoff,
   - the edge weight is the estimated residue-level contact area;

9. Identifying hydrophobic clusters as weakly connected components in the residue
   contact graph, so that residues connected directly or indirectly by sufficient
   hydrophobic contact are grouped into the same cluster;

10. Calculating the total area of each hydrophobic cluster as the sum of all
    residue-to-residue contact edge areas inside that cluster;

11. Generating a structured JSON report containing all detected hydrophobic
    clusters and their member residues.
'''

# dependency:
'''
Biopython
NumPy
SciPy
NetworkX
'''

# reference:
'''
The hydrophobic cluster calculation is adapted from the Protlego implementation:
https://github.com/Hoecker-Lab/protlego/blob/master/protlego/structural/clusters.py

Original Protlego reference:
Flores SC, et al. Protlego: A Python package for the analysis and design of
chimeric proteins.
https://github.com/Hoecker-Lab/protlego
'''