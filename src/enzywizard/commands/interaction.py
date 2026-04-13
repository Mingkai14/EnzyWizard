from __future__ import annotations

from argparse import Namespace
from ..services.interaction_service import run_interaction_service


def add_interaction_parser(subparsers) -> None:
    parser = subparsers.add_parser("interaction",help="Calculate protein-substrate and intra-protein interactions for an input CIF/PDB structure using docked substrate SDF files.")
    parser.add_argument("-i","--input_path",required=True,help="Path to input CIF/PDB file.")
    parser.add_argument("-s","--substrate_names",required=False,default=None,help="Input substrate names separated by ','. Each substrate name must match the corresponding docked SDF file name in substrate_dir. If omitted together with --substrate_dir, only intra-protein interactions will be calculated.")
    parser.add_argument("-d","--substrate_dir",required=False,default=None,help="Optional path to a directory containing docked substrate SDF files. Must be provided together with --substrate_names.")
    parser.add_argument("-o","--output_dir",required=True,help="Path to a directory for outputting the interaction JSON report.")
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


# input parameters:
'''
-i --input_path Required. Path to input cleaned protein structure file (CIF or PDB).

-s --substrate_names Required. Input substrate names separated by ','.

Examples:
- docked_glucose
- docked_glucose,docked_fructose
- ...

Each substrate name must match the corresponding docked substrate SDF file name
in substrate_dir.

This parameter represents the single docked substrate/multiple simultaneously docked substrates used for interaction analysis.

For input substrate_names, the program searches
substrate_dir for matched docked substrate SDF files and loads them as
substrate Mol(3D) objects for interaction calculation.

-d --substrate_dir Required. Path to a directory containing docked substrate
SDF files.

The program only reads matched .sdf files from this directory.

-o --output_dir Required. Output directory for saving the interaction JSON report.

--hbond_da_max_distance Optional. Maximum donor-acceptor distance cutoff for
hydrogen bond detection (default: 3.9).

--hbond_ha_max_distance Optional. Maximum hydrogen-acceptor distance cutoff for
hydrogen bond detection (default: 2.5).

--hbond_angle Optional. Minimum donor-hydrogen-acceptor angle cutoff for
hydrogen bond detection in degrees (default: 90.0).

--ionic_distance_cutoff Optional. Maximum distance cutoff for ionic interaction
detection (default: 4.0).

--ppstack_center_distance_cutoff Optional. Maximum aromatic ring-center
distance cutoff for pi-pi stacking detection (default: 6.5).

--pication_distance_cutoff Optional. Maximum ring-cation distance cutoff for
pi-cation interaction detection (default: 5.0).

--pication_angle_cutoff Optional. Maximum angle cutoff for pi-cation
interaction detection in degrees (default: 45.0).

--ssbond_max_distance Optional. Maximum sulfur-sulfur distance cutoff for
disulfide bond detection (default: 2.5).
'''

# output content:
'''
The program outputs:

1. A JSON report containing:
   - "output_type": "enzywizard_interaction"
   - "interactions": all detected interaction records
   - "interactions_statistics": summarized interaction counts

1.1 Each item in "interactions" contains:
   - "interaction": interaction type
   - "node1": node 1 information
   - "node2": node 2 information

   Supported interaction types include:
   - "HBOND"
   - "IONIC"
   - "VDW"
   - "PIPISTACK"
   - "PICATION"
   - "SSBOND"

1.2 Amino-acid node format:
   - "aa_index"
   - "aa_name"
   - "node_type": "amino_acid"

1.3 Substrate node format:
   - "substrate_index"
   - "substrate_name"
   - "node_type": "substrate"

1.4 "interactions_statistics" contains three scopes:
   - "overall"
   - "intra_protein"
   - "protein_substrate"

For each scope, the report contains:
   - "count": total number of detected interactions for each interaction type
   - "unique_pair_count": number of unique node pairs for each interaction type

'''

# functionality/process
'''
It performs interaction analysis by:

1. Reading the cleaned input protein structure from input_path;

2. Checking whether the structure satisfies cleaned-structure requirements;

3. Loading the same input structure as an OpenMM Modeller object for
   coordinate-based interaction calculation;

4. Checking whether the protein structure contains hydrogen atoms;

5. Parsing substrate_names using ',' as the separator to obtain one or more
   substrate names for interaction analysis;

6. Searching substrate_dir for matched docked substrate SDF files and loading
   them as substrate Mol(3D) objects with explicit hydrogen checking;

7. Filtering invalid substrate molecules and removing substrate molecules that
   are not spatially docked to the protein;

8. Detecting six types of interactions:
    - hydrogen bond (HBOND)
    - ionic interaction (IONIC)
    - van der Waals contact (VDW)
    - pi-pi stacking (PIPISTACK)
    - pi-cation interaction (PICATION)
    - disulfide bond (SSBOND)

9. Merging and sorting all detected interactions into a interaction list with
    standardized node information;

10. Summarizing interaction counts for:
    - all interactions
    - intra-protein interactions
    - protein-substrate interactions

11. Writing a structured JSON report.
'''

# algorithm details
'''
The program detects interactions using the following algorithm design.

1. Common preprocessing

1.1 Protein coordinates are read from OpenMM Modeller.

1.2 OpenMM positions are converted from nm to Å.

1.3 A substrate is considered docked when at least one protein heavy atom and
   one substrate heavy atom are within 6.5 Å.

1.4 For intra-protein interactions, residue pairs with too small 
sequence separation are filtered using min_residue_index_gap = 3 
to reduce local contacts.



2 Hydrogen bond detection (HBOND)

2.1 Protein donor and acceptor atoms are defined from:
   - backbone atoms
   - side-chain functional groups

2.2 Substrate donor and acceptor atoms are identified using RDKit
   chemical features.

2.3 A hydrogen bond is accepted when all conditions are satisfied:
   - D-A distance <= hbond_da_max_distance (default: 3.9 Å)
   - H-A distance <= hbond_ha_max_distance (default: 2.5 Å)
   - angle(D-H-A) >= hbond_angle (default: 90°)

2.4 Intra-protein:
   - same residue and local contacts are excluded

2.5 Protein-substrate:
   - donor/acceptor roles are considered in both directions
   - multiple detections between the same residue and substrate are merged


3 Ionic interaction detection (IONIC)

3.1 Protein charged centers are defined as:
   - ASP/GLU → carboxylate center (anion)
   - LYS/ARG/HIS → side-chain center (cation)

3.2 Substrate charged centers are defined using RDKit formal charges.

3.3 Oppositely charged center pairs are evaluated:
   - protein–protein
   - protein–substrate

3.4 An ionic interaction is accepted when:
   - center distance <= ionic_distance_cutoff (default: 4.0 Å)

3.5 Intra-protein:
   - same residue and local contacts are excluded

3.6 Protein-substrate:
   - repeated detections are merged

3.7 Output interaction type:
   - "IONIC"


4 van der Waals contact detection (VDW)


4.1 Protein and substrate atoms are assigned van der Waals radii
   based on element type.

4.2 Atom pairs are evaluated using distance vs radius sum:

   - let r1 + r2 = sum of radii
   - contact if:
     (1 - mu)*(r1 + r2) <= d <= (1 + mu)*(r1 + r2)

   - default: mu = 0.01

4.3 KD-tree is used to accelerate neighbor search.

4.4 Intra-protein:
   - same residue and local contacts are excluded

4.5 Protein-substrate:
   - repeated detections are merged


5 pi-pi stacking detection (PIPISTACK)


5.1 Protein aromatic rings are defined from residue templates:
   - PHE, TYR, HIS, TRP

5.2 Substrate aromatic rings are detected using RDKit:
   - ring size >= 5
   - all atoms aromatic

5.3 For each ring:
   - center = average coordinate
   - normal = fitted plane normal

5.4 Ring pairs are accepted when:
   - center distance <= ppstack_center_distance_cutoff (default: 6.5 Å)

5.5 Additional angular features are computed for geometry classification
   (parallel, T-shaped, etc.), but not used for filtering.

5.6 Intra-protein:
   - same residue and local contacts are excluded

5.7 Protein-substrate:
   - repeated detections are merged


6 pi-cation interaction detection (PICATION)


6.1 Protein aromatic rings are defined as in pi-pi stacking.

6.2 Protein cation centers:
   - LYS (NZ)
   - ARG (guanidinium center)

6.3 Substrate cation centers are defined using RDKit formal charges.

6.4 Ring–cation pairs are accepted when:
   - distance <= pication_distance_cutoff (default: 5.0 Å)
   - angle between ring normal and ring–cation vector <= pication_angle_cutoff
     (default: 45°)

6.5 Intra-protein:
   - same residue and local contacts are excluded

6.6 Protein-substrate:
   - repeated detections are merged


7 Disulfide bond detection (SSBOND)


7.1 Cysteine SG atoms are collected from protein.

7.2 SG–SG pairs are evaluated:
   - distance <= ssbond_max_distance (default: 2.5 Å)

7.3 Intra-protein:
   - same residue and local contacts are excluded

'''

# dependency:
'''
RDKit
OpenMM
Biopython
'''


# reference:
'''
Martin et al., RING 4.0: residue interaction network generation for protein structures and ensembles
https://ring.biocomputingup.it/about

Piovesan et al., RING 2.0: fast generation of residue interaction networks
https://doi.org/10.1093/bioinformatics/btw203

Bondi, A. van der Waals Volumes and Radii
J. Phys. Chem. 1964, 68, 3, 441–451
https://doi.org/10.1021/j100785a001

RDKit: Open-source cheminformatics
https://www.rdkit.org/

OpenMM: Molecular simulation toolkit
https://openmm.org/

Biopython: Structural bioinformatics tools
https://biopython.org/
'''