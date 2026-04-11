from __future__ import annotations
from argparse import Namespace
from ..services.dock_service import run_dock_service


def add_dock_parser(subparsers) -> None:
    parser = subparsers.add_parser("dock",help="Perform substrate docking for an input CIF/PDB structure using substrate SDF files.")
    parser.add_argument("-i", "--input_path",required=True,help="Path to input CIF/PDB file.")
    parser.add_argument("-s", "--substrate_names",required=True,help="Input substrate names separated by ','. Each substrate name must match corresponding SDF file names in substrate_dir.")
    parser.add_argument("-d","--substrate_dir",required=True,help="Path to a directory containing input substrate SDF files used for docking.")
    parser.add_argument("-o", "--output_dir",required=True,help="Path to a directory for outputting docked substrate SDF files, docked protein-substrate complex CIF/PDB files, and a JSON report.")
    parser.add_argument("--max_docking_attempt_num",type=int,default=20,help="Maximum number of docking attempts (default: 20).")
    parser.add_argument("--early_stop", type=lambda x: str(x).lower() in ["true", "1", "yes"], default=False, help="Whether to stop immediately after the first successful docking result (True/False, default: False). If True, the returned result may not be the global best.")
    parser.add_argument("--exhaustiveness",type=int,default=16,help="Exhaustiveness of AutoDock Vina search (default: 16). Larger values may improve docking search coverage but increase runtime.")
    parser.add_argument("--cpu",type=int,default=0,help="Number of CPUs used by AutoDock Vina (default: 0). A value of 0 lets Vina decide automatically.")
    parser.add_argument("--min_rad",type=float,default=1.8,help="Minimum probe radius used in pocket detection (default: 1.8). Smaller values may detect narrower cavities, but overly small values may cause PyVOL/MSMS failure.")
    parser.add_argument("--max_rad",type=float,default=6.2,help="Maximum probe radius used in pocket detection (default: 6.2). Larger values may detect broader cavities, but overly large values may cause PyVOL/MSMS failure.")
    parser.add_argument("--min_volume",type=int,default=50,help="Minimum pocket volume threshold used to filter detected pocket regions (default: 50). Larger values retain only larger pocket candidates.")

    parser.set_defaults(func=run_dock)


def run_dock(args: Namespace) -> None:
    run_dock_service(
        input_path=args.input_path,
        substrate_names=args.substrate_names,
        substrate_dir=args.substrate_dir,
        output_dir=args.output_dir,
        max_docking_attempt_num=args.max_docking_attempt_num,
        early_stop=args.early_stop,
        exhaustiveness=args.exhaustiveness,
        cpu=args.cpu,
        min_rad=args.min_rad,
        max_rad=args.max_rad,
        min_volume=args.min_volume,
    )


# input parameters:
'''
-i --input_path Required. Path to input cleaned protein structure file (CIF or PDB).

-s --substrate_names Required. Input substrate names separated by ','.

Examples:
- glucose
- glucose,fructose
- acetate,ethanol

This parameter represents a multi-substrate combination for docking.

This parameter helps program to support:
- single-substrate docking
- multi-substrate simultaneous docking
- multi-conformation substrate combining

For intput substrate_names "SubstrateA,SubstrateB", the program searches substrate_dir for matched SDF files:

- SubstrateA.sdf
- SubstrateA_1.sdf
- SubstrateA_2.sdf
- SubstrateA_3.sdf
- ...

- SubstrateB.sdf
- SubstrateB_1.sdf
- SubstrateB_2.sdf
- SubstrateB_3.sdf
- ...

The matched SDF files are treated as different conformations to combine and attempt docking:

- SubstrateA_1.sdf + SubstrateB_1.sdf
- SubstrateA_1.sdf + SubstrateB_2.sdf
- SubstrateA_1.sdf + SubstrateB_3.sdf
- ...

Duplicate substrate names are not allowed.

-d --substrate_dir Required. Path to a directory containing input substrate SDF files. The program only reads matched .sdf files from this directory.

-o --output_dir Required. Output directory for saving docked substrate files, docked protein-substrate complex files, and a JSON report.

--max_docking_attempt_num Optional. Maximum number of docking attempts (default: 20).

--early_stop Optional. Whether to stop immediately after the first successful docking result (True/False, default: False).

If True, the program returns the first successful result and does not continue searching for a better one.

--exhaustiveness Optional. Exhaustiveness of AutoDock Vina search (default: 16). Larger values may improve search coverage but increase runtime.

--cpu Optional. Number of CPUs used by AutoDock Vina (default: 0). A value of 0 lets Vina decide automatically.

--min_rad Optional. Minimum probe radius used for pocket detection (default: 1.8). Smaller values may detect narrower cavities, but overly small values may cause PyVOL/MSMS failure.

--max_rad Optional. Maximum probe radius used for pocket detection (default: 6.2). Larger values may detect broader cavities, but overly large values may cause PyVOL/MSMS failure.

--min_volume Optional. Minimum pocket volume threshold used to filter detected pockets (default: 50). Larger values retain only larger pocket candidates.
'''


# output content:
'''
The program outputs:

1. A JSON report containing:
   - "output_type": "enzywizard_dock"
   - "docked_result": the best docking result record

1.1 The "docking_result" includes:
   - complex_name
   - docking_score
   - substrate_names
   - docking_box_center
   - docking_box_size
   - docked_substrates

1.2 The "docked_substrates" includes a list of docked substrate records, and each record contains:
   - substrate_name
   - conformation_name
   - docked_center_coord

2. Files:
   - one docked SDF file for each docked substrate
   - one docked protein-substrate complex CIF file
   - one JSON report

'''


# functionality/process
'''
It performs docking by:

1. Reading the cleaned input protein structure from input_path;

2. Checking whether the structure satisfies cleaned-structure requirements;

3. Detecting pocket regions from the protein structure using PyVOL;

4. Calculating a global docking box from the whole protein structure;

5. Parsing substrate_names using ',' as the separator to obtain one or more
   substrate names for docking;

6. Searching substrate_dir for matched SDF files for each substrate name from input "substrate_names" (e.g. "SubstrateA,SubstrateB"):

- SubstrateA.sdf
- SubstrateA_1.sdf
- SubstrateA_2.sdf
- SubstrateA_3.sdf
- ...

- SubstrateB.sdf
- SubstrateB_1.sdf
- SubstrateB_2.sdf
- SubstrateB_3.sdf
- ...

7. Treating matched SDF files of the same substrate as alternative
   conformations, and automatically enumerating substrate
   combinations for docking;

8. Converting the protein structure into receptor PDBQT format and converting
   each matched substrate SDF into ligand PDBQT format;

9. Building docking boxes from:
   - detected pocket boxes
   - one global whole-structure box

10. Iterating over substrate combinations and docking boxes, and performing
    AutoDock Vina docking for single substrate docking/simultaneous multi-substrate docking;

11. Reading docking poses and energies returned by Vina;

12. Writing docked substrate SDF files;

13. Writing the docked protein-substrate complex CIF file;

14. Selecting and saving the best docking result, then generating a structured JSON report.
'''


# dependency:
'''
AutoDock Vina
Meeko (for PDBQT preparation)
RDKit
Biopython
PyVOL
MSMS (used by PyVOL)
'''


# reference:
'''
Eberhardt et al., AutoDock Vina 1.2.0: New docking methods, expanded force field, and Python bindings
https://doi.org/10.1021/acs.jcim.1c00203

Trott & Olson, AutoDock Vina: improving the speed and accuracy of docking
https://doi.org/10.1002/jcc.21334

AutoDock Vina official documentation
https://vina.scripps.edu/

Meeko: Preparation of small molecules for AutoDock
https://github.com/forlilab/Meeko

RDKit: Open-source cheminformatics
https://www.rdkit.org/

Biopython: Structural bioinformatics tools
https://biopython.org/

PyVOL: Python interface for binding pocket detection
https://github.com/schlessinger-lab/pyvol

MSMS: Molecular surface calculation
Sanner et al., Reduced surface: an efficient way to compute molecular surfaces
https://doi.org/10.1002/jcc.540150805
'''