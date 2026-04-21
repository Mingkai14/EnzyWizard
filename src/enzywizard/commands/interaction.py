from __future__ import annotations

from argparse import Namespace
from ..services.interaction_service import run_interaction_service


def add_interaction_parser(subparsers) -> None:
    parser = subparsers.add_parser("interaction",help="Calculate protein/protein-substrate interactions and generate a detailed JSON report.")
    parser.add_argument("-i","--input_path",required=True,help="Path to the input cleaned protein structure file in CIF or PDB format.")
    parser.add_argument("-s","--substrate_names",required=False,default=None,help="Input substrate names separated by ','. Each substrate name must match the corresponding docked SDF file name in substrate_dir. If omitted together with --substrate_dir, only intra-protein interactions will be calculated.")
    parser.add_argument("-d","--substrate_dir",required=False,default=None,help="Optional path to a directory containing docked substrate SDF files. Must be provided together with --substrate_names.")
    parser.add_argument("-o","--output_dir",required=True,help="Directory to save the interaction JSON report.")
    parser.add_argument("--hbond_da_max_distance",type=float,default=3.9,help="Maximum donor-acceptor distance cutoff for hydrogen bond detection (default: 3.9).")
    parser.add_argument("--hbond_ha_max_distance",type=float,default=2.5,help="Maximum hydrogen-acceptor distance cutoff for hydrogen bond detection (default: 2.5).")
    parser.add_argument("--hbond_angle",type=float,default=90.0,help="Minimum donor-hydrogen-acceptor angle cutoff for hydrogen bond detection (default: 90.0).")
    parser.add_argument("--ionic_distance_cutoff",type=float,default=4.0,help="Maximum distance cutoff for ionic bond detection (default: 4.0).")
    parser.add_argument("--ppstack_center_distance_cutoff",type=float,default=6.5,help="Maximum ring-center distance cutoff for pi-pi stacking detection (default: 6.5).")
    parser.add_argument("--pication_distance_cutoff",type=float,default=5.0,help="Maximum ring-cation distance cutoff for pi-cation interaction detection (default: 5.0).")
    parser.add_argument("--pication_angle_cutoff",type=float,default=45.0,help="Maximum angle cutoff for pi-cation interaction detection (default: 45.0).")
    parser.add_argument("--ssbond_max_distance",type=float,default=2.5,help="Maximum sulfur-sulfur distance cutoff for disulfide bond detection (default: 2.5).")

    parser.set_defaults(func=run_interaction)


def run_interaction(args: Namespace) -> None:
    run_interaction_service(
        input_path=args.input_path,
        output_dir=args.output_dir,
        substrate_names=args.substrate_names,
        substrate_dir=args.substrate_dir,
        da_max_distance_A=args.hbond_da_max_distance,
        ha_max_distance_A=args.hbond_ha_max_distance,
        dha_min_angle_deg=args.hbond_angle,
        ionic_distance_cutoff_A=args.ionic_distance_cutoff,
        ring_center_distance_cutoff_A=args.ppstack_center_distance_cutoff,
        ring_cation_distance_cutoff_A=args.pication_distance_cutoff,
        ring_cation_angle_cutoff_deg=args.pication_angle_cutoff,
        ss_max_distance_A=args.ssbond_max_distance,
    )


# ==============================
# Command: enzywizard-interaction
# ==============================

# brief introduction:
'''
EnzyWizard-Interaction is a command-line tool for calculating protein/protein-substrate
interactions and generating a detailed JSON report.
It takes a cleaned CIF or PDB protein structure and a substrate directory as input and supports 
protein-(multi)substrate(s) interaction analysis. The tool detects six types of interactions:
hydrogen bond (HBOND), ionic interaction (IONIC), van der Waals contact (VDW), pi-pi stacking (PIPISTACK),
pi-cation interaction (PICATION), disulfide bond (SSBOND).
The geometric definitions and interaction calculation formula follow
the residue interaction network methods RING4.
To improve efficiency, the program uses vectorized geometry operations and
KD-tree accelerated neighbour searching where applicable, enabling fast
interaction detection on protein / docked protein-substrate systems.
The tool outputs structured interaction records and summarized statistics
suitable for downstream applications such as enzyme-substrate interaction
analysis and residue interaction graph construction.
'''

# example usage:
'''
Example command:

enzywizard-interaction -i examples/input/cleaned_3GP6.cif -s docked_glucose,docked_fructose -d examples/input/ -o examples/output/
'''

# input parameters:
'''
-i, --input_path
Required.
Path to the input cleaned protein structure file in CIF or PDB format.

-s, --substrate_names
Optional.
Input substrate names separated by ','.

Examples:
- docked_glucose
- docked_glucose,docked_fructose

Each substrate name must match the corresponding docked substrate SDF file name
in substrate_dir.

If --substrate_names and --substrate_dir are both omitted, the program only
calculates intra-protein interactions.

-d, --substrate_dir
Optional.
Path to a directory containing docked substrate SDF files.

This parameter must be provided together with --substrate_names.

The program only reads matched .sdf files from this directory.

-o, --output_dir
Required.
Directory to save the interaction JSON report.

--hbond_da_max_distance
Optional.
Maximum donor-acceptor distance cutoff for hydrogen bond detection.

Default:
  3.9

--hbond_ha_max_distance
Optional.
Maximum hydrogen-acceptor distance cutoff for hydrogen bond detection.

Default:
  2.5

--hbond_angle
Optional.
Minimum donor-hydrogen-acceptor angle cutoff for hydrogen bond detection in degrees.

Default:
  90.0

--ionic_distance_cutoff
Optional.
Maximum distance cutoff for ionic interaction detection.

Default:
  4.0

--ppstack_center_distance_cutoff
Optional.
Maximum aromatic ring-center distance cutoff for pi-pi stacking detection.

Default:
  6.5

--pication_distance_cutoff
Optional.
Maximum ring-cation distance cutoff for pi-cation interaction detection.

Default:
  5.0

--pication_angle_cutoff
Optional.
Maximum angle cutoff for pi-cation interaction detection in degrees.

Default:
  45.0

--ssbond_max_distance
Optional.
Maximum sulfur-sulfur distance cutoff for disulfide bond detection.

Default:
  2.5
'''

# output content:
'''
The program outputs the following file into the output directory:

1. A JSON report
   - interaction_report_{protein_name}.json

   or

   - interaction_report_{protein_name}_{substrate_names}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_interaction"

   - "interactions"
     A list describing all detected interaction records.

     Each entry contains:
     - "interaction"
       Interaction type.

       Supported interaction types include:
       - "HBOND"
       - "IONIC"
       - "VDW"
       - "PIPISTACK"
       - "PICATION"
       - "SSBOND"

     - "node1"
       Information of the first interaction node.

     - "node2"
       Information of the second interaction node.

     Node information is stored in one of the following formats.

     Amino-acid node:
     - "aa_index"
       Residue index in the cleaned structure.

     - "aa_name"
       Residue one-letter amino acid code.

     - "node_type"
       "amino_acid"

     Substrate node:
     - "substrate_index"
       Internal substrate index in the current interaction calculation.

     - "substrate_name"
       Input substrate name.

     - "node_type"
       "substrate"

   - "interactions_statistics"
     A dictionary summarizing detected interaction counts.

     It contains three scopes:
     - "overall"
     - "intra_protein"
     - "protein_substrate"

     For each scope, the report contains:
     - "count"
       Total number of detected interactions for each interaction type.

     - "unique_pair_count"
       Number of unique node pairs for each interaction type.
'''

# Process:
'''
This command processes the input cleaned protein structure as follows:

1. Load the input structure
   - Read the cleaned CIF or PDB file using Biopython (Bio.PDB).
   - Resolve the protein name from the input filename.

2. Validate input conditions
   - Check that the input file exists.
   - Validate that the structure satisfies the cleaned-structure requirement.
   - Validate interaction cutoff parameters.

3. Load OpenMM modeller
   - Read the same input structure as an OpenMM Modeller object.
   - Build coordinate-based protein representation for geometric interaction calculation.

4. Check hydrogen completeness
   - Detect whether the protein structure contains hydrogen atoms.
   - Report warnings if hydrogen atoms are missing or too few.

5. Parse substrate inputs
   - Determine whether the program runs in:
     - protein-only mode
     - protein-substrate mode

6. Load docked substrates
   - If substrate_names and substrate_dir are provided, search substrate_dir for matched docked substrate SDF files.
   - Load matched SDF files as substrate Mol(3D) objects.

7. Filter substrate molecules
   - Check whether each substrate molecule is a valid 3D structure.
   - Check whether explicit hydrogen atoms are present.
   - Check whether the substrate is spatially docked to the protein.
   - Retain valid substrates for interaction calculation.

8. Build geometric tables
   - Convert OpenMM coordinates from nm to Å.
   - Build protein atom/bond tables.
   - Build substrate atom/bond tables.

9. Detect hydrogen bonds
   - Identify protein donor and acceptor atoms from backbone and side-chain functional groups.
   - Identify substrate donor and acceptor atoms using RDKit chemical features.
   - Detect hydrogen bonds using distance and angle constraints.

10. Detect ionic interactions
   - Build charged centers for protein residues and substrate atoms.
   - Detect oppositely charged center pairs within the distance cutoff.

11. Detect van der Waals contacts
   - Assign van der Waals radii based on atom element type.
   - Detect atom pairs whose distances satisfy radius-sum contact constraints.

12. Detect pi-pi stacking
   - Identify protein aromatic rings from residue templates.
   - Identify substrate aromatic rings using RDKit aromatic ring information.
   - Detect ring pairs using ring-center distance cutoff.

13. Detect pi-cation interactions
   - Identify aromatic rings and cation centers from protein and substrate.
   - Detect ring-cation pairs using distance and angle constraints.

14. Detect disulfide bonds
   - Identify cysteine SG atoms from protein.
   - Detect sulfur-sulfur pairs using distance cutoff.

15. Merge interaction records
   - Convert all detected raw edges into standardized interaction records.
   - Build node information for amino-acid nodes and substrate nodes.
   - Sort interaction records into a stable output order.

16. Summarize interaction statistics
   - Count detected interactions for:
     - all interactions
     - intra-protein interactions
     - protein-substrate interactions
   - Calculate unique node-pair counts for each interaction type.

17. Save outputs
   - Generate and save a JSON report containing interaction records and summary statistics.
'''

# dependencies:
'''
- RDKit
- OpenMM
- Biopython
- NumPy
- SciPy
'''

# references:
'''
- Martin et al., RING 4.0: residue interaction network generation for protein structures and ensembles
  https://ring.biocomputingup.it/about

- Piovesan et al., RING 2.0: fast generation of residue interaction networks
  https://doi.org/10.1093/bioinformatics/btw203

- Bondi, A. van der Waals Volumes and Radii
  J. Phys. Chem. 1964, 68, 3, 441–451
  https://doi.org/10.1021/j100785a001

- RDKit:
  https://www.rdkit.org/

- OpenMM:
  https://openmm.org/

- Biopython:
  https://biopython.org/
'''