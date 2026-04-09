from __future__ import annotations
from argparse import Namespace
from ..services.dock_service import run_dock_service


def add_dock_parser(subparsers) -> None:
    parser = subparsers.add_parser("dock",help="Perform substrate docking for an input CIF/PDB structure using substrate SDF files.")
    parser.add_argument("-i", "--input_path",required=True,help="Path to input CIF/PDB file.")
    parser.add_argument("-s", "--substrate_names",required=True,help="Input substrate names separated by ','. Each substrate name must match corresponding SDF file names in substrate_dir.")
    parser.add_argument("--substrate_dir",required=True,help="Path to a directory containing input substrate SDF files used for docking.")
    parser.add_argument("-o", "--output_dir",required=True,help="Path to a directory for outputting docked substrate SDF files, docked protein-substrate complex CIF/PDB files, and a JSON report.")
    parser.add_argument("--max_docking_result_num",type=int,default=3,help="Maximum number of successful docking results to keep (default: 3). Once this number is reached, the workflow stops and returns the currently collected best-scoring results.")
    parser.add_argument("--max_docking_attempt_num",type=int,default=20,help="Maximum number of docking attempts (default: 20).")
    parser.add_argument("--max_pose_read_num",type=int,default=1,help="Maximum number of poses read from each successful Vina docking run (default: 1). Larger values allow multiple top-ranked poses to be retained from a single docking attempt.")
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
        max_docking_result_num=args.max_docking_result_num,
        max_docking_attempt_num=args.max_docking_attempt_num,
        max_pose_read_num=args.max_pose_read_num,
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
- substrate1,substrate2,substrate3

This parameter supports both single-substrate docking and multi-substrate
simultaneous docking.

Each substrate name must match corresponding SDF file names in substrate_dir.

For a substrate name such as glucose, the program searches substrate_dir for:
   - glucose.sdf
   - glucose_1.sdf
   - glucose_2.sdf
   - glucose_3.sdf
   - ...

If multiple matched SDF files exist for one substrate, they are treated as
alternative input conformations/protomers of the same substrate.

For multi-substrate docking, the program selects one matched SDF file from each
candidate conformations group of one substrate at a time and combines them into one multi-substrate docking input set.

For example, if:
- glucose has glucose_1.sdf and glucose_2.sdf
- fructose has fructose_1.sdf and fructose_2.sdf

the program will try multi-substrate docking for combinations such as:
- glucose_1 + fructose_1
- glucose_1 + fructose_2
- glucose_2 + fructose_1
- glucose_2 + fructose_2

Duplicate substrate names are not allowed.

--substrate_dir Required. Path to a directory containing input substrate SDF files.

The program reads substrate SDF files only from this directory.

For substrate_name = glucose:
accepted file names include:
  glucose.sdf
  glucose_1.sdf
  glucose_2.sdf
  ...

- file extensions must be .sdf
- unrelated SDF files are ignored
- each matched SDF file is treated as one candidate input conformation for docking

-o --output_dir Required. Output directory for saving docked structure files
and the JSON report.

--max_docking_result_num Optional. Maximum number of successful docking
results to keep (default: 3).

This parameter controls the early stopping threshold of the docking workflow.
Once enough successful docking results have been collected, the workflow stops
and returns the currently collected top-scoring results.

--max_docking_attempt_num Optional. Maximum number of docking attempts
(default: 20).

Each Vina docking call for one substrate-combination under one docking box
counts as one attempt. This parameter limits total search cost.

--max_pose_read_num Optional. Maximum number of poses read from each
successful docking run (default: 1).

This parameter controls how many top-ranked poses are extracted from one Vina
run. Larger values may retain more candidate binding modes per docking call.

--exhaustiveness Optional. Exhaustiveness of AutoDock Vina search
(default: 16).

This parameter controls the search depth of Vina. Larger values may improve
search completeness but also increase runtime.

--cpu Optional. Number of CPUs used by AutoDock Vina (default: 0).

A value of 0 allows Vina to determine CPU usage automatically.

--min_rad Optional. Minimum probe radius used for pocket detection
(default: 1.8).

This parameter controls the smallest probe sphere used during pocket/cavity
detection. Smaller values may detect narrower cavities, but overly small values
may cause PyVOL/MSMS errors.

--max_rad Optional. Maximum probe radius used for pocket detection
(default: 6.2).

This parameter controls the largest probe sphere used during pocket/cavity
expansion. Larger values may detect broader pockets, but overly large values
may cause PyVOL/MSMS errors.

--min_volume Optional. Minimum pocket volume threshold used to filter detected
pockets (default: 50).

This parameter controls the lower size limit of retained pocket candidates.
Larger values remove smaller pocket regions from docking box generation.
'''


# output content:
'''
The program outputs:

1. A JSON report recording:
   - "output_type": "enzywizard_dock"
   - "docked_results": a list of docking result records

2. For each docking result record, the JSON report includes:
   - complex_name
   - docking_score
   - substrate_names
   - docking_box_center
   - docking_box_size
   - pose_index
   - docked_substrates

3. "docked_substrates" is a list of docked substrate records, where each
   entry contains:
   - substrate_name
   - docked_center_coord
   - docked_sdf_path

4. Each retained docking result corresponds to one successful simultaneous
   docking output for one substrate combination under one docking box and
   one docking pose.

5. For each retained docking result, the program writes:
   - one docked SDF file for each docked substrate
   - one docked protein-substrate complex CIF file
   - one docked protein-substrate complex PDB file

6. For multi-substrate docking, one result contains all docked substrates
   belonging to the same docking run.

7. Docked output file naming follows the internal docking result index, for
   example:
   - docked1_glucose.sdf
   - docked1_fructose.sdf
   - docked1_protein_glucose_fructose.cif
   - docked1_protein_glucose_fructose.pdb
   - docked2_glucose.sdf
   - ...

'''


# functionality/process
'''
It performs docking by:

1. Reading the cleaned input protein structure from input_path;

2. Checking whether the structure satisfies cleaned-structure requirements
   before docking;

3. Detecting pocket regions from the protein structure using PyVOL;

4. Calculating a global structure box from the whole protein structure;

5. Parsing substrate_names using ',' as the separator to obtain one or more
   substrate names for docking;

6. Searching substrate_dir for matching SDF files for each substrate name:
   - substrate.sdf
   - substrate_1.sdf
   - substrate_2.sdf
   - ...
   These files are treated as alternative candidate conformations/protomers
   for the same substrate;

7. Grouping matched SDF files by substrate name, then constructing docking
   combinations by selecting one candidate SDF file from each substrate group;

8. Converting the protein structure into receptor PDBQT format using Meeko;
   converting each substrate SDF file into ligand PDBQT format using Meeko;

9. Building docking boxes from:
    - detected pocket boxes
    - one global whole-structure box

10. Iterating over substrate combinations and docking boxes, and performing
    AutoDock Vina docking for simultaneous multi-substrate docking;

11. Reading top docking poses and energies returned by Vina;

12. Reconstructing docked substrate 3D molecules and saving them as SDF files;

13. Writing docked protein-substrate complex structures as CIF and PDB files;

14. Collecting all successful docking results, sorting them by docking score,
    and generating a structured JSON report.
'''


# algorithm:
'''
Core algorithms and methods used in this workflow include:

1. Pocket detection:
   - Pocket candidates are computed by PyVOL.
   - Pocket boxes are derived from detected pocket center coordinates and
     box boundaries.
   - A global whole-protein box is also generated as an additional docking box.

2. Multi-substrate and multi-conformation combination strategy:
   - For each requested substrate, one matched SDF file is selected at a time.
   - Multiple matched SDF files for the same substrate are treated as
     alternative conformations/protomers of that substrate.
   - The program enumerates combinations by taking one candidate SDF from
     each substrate group.
   - Each combination represents one candidate multi-substrate input set
     for simultaneous docking.

3. Ligand/receptor preparation:
   - Protein structure is converted to receptor PDBQT by Meeko.
   - Substrate SDF structures are converted to ligand PDBQT by Meeko.
   - Ligand atom index mapping is retained using Meeko's index map output.

4. Docking engine:
   - AutoDock Vina is used for docking calculation.
   - Vina maps are computed for each docking box.
   - Docking is run with user-controlled exhaustiveness, CPU count, and number
     of poses to retain.
   - Multiple ligands in one substrate combination are docked simultaneously
     in the same Vina run.

5. Pose reconstruction:
   - Vina output pose strings are split into MODEL blocks.
   - Ligand blocks are extracted from each pose.
   - Docked atom coordinates are recovered from pose ligand blocks.

6. Docked structure writing:
   - Recovered docked coordinates are written back into RDKit molecule objects.
   - Docked ligands are saved as SDF files.
   - Protein and docked ligands are merged into one complex structure and saved
     as CIF/PDB files.

7. Result ranking:
   - Successful docking results are ranked by Vina docking score
     (lower score is better rank).
   - The workflow stops early if max_docking_result_num is reached, or stops
     when max_docking_attempt_num is exceeded.
'''