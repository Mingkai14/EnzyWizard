from __future__ import annotations
from argparse import Namespace
from ..services.pocket_service import run_pocket_service


def add_pocket_parser(subparsers) -> None:
    parser = subparsers.add_parser("pocket",help="Detect and characterize binding pockets from a cleaned protein structure and generate a detailed JSON report.")
    parser.add_argument("-i", "--input_path", required=True,help="Path to the input cleaned protein structure file in CIF or PDB format.")
    parser.add_argument("-o", "--output_dir",required=True,help="Directory to save the JSON report.")
    parser.add_argument("--min_rad",type=float,default=1.8,help="Minimum probe radius used by PyVOL for cavity detection (default: 1.8). Smaller values allow detection of narrower cavities, but excessively small values may lead to PyVOL failure.")
    parser.add_argument("--max_rad",type=float,default=6.2,help="Maximum probe radius used by PyVOL for cavity detection (default: 6.2). Larger values allow detection of broader pockets, but excessively large values may lead to PyVOL failure.")
    parser.add_argument("--min_volume",type=int,default=50,help="Minimum pocket volume threshold (default: 50). Pockets with volume below this value will be discarded.")

    parser.set_defaults(func=run_pocket)


def run_pocket(args: Namespace) -> None:
    run_pocket_service(input_path=args.input_path,output_dir=args.output_dir,min_rad=args.min_rad,max_rad=args.max_rad,min_volume=args.min_volume)

# ==============================
# Command: enzywizard-pocket
# ==============================

# brief introduction:
'''
EnzyWizard-Pocket is a command-line tool for detecting and characterizing
binding pockets from a cleaned protein structure and generating a detailed JSON report.
It takes a CIF or PDB file as input and identifies potential cavity regions using 
PyVOL, a geometry-based pocket detection program. PyVOL detects pockets by 
simulating probe spheres rolling over the protein surface,
identifying cavities based on geometric accessibility and spatial continuity.
Each detected pocket is represented as a cluster of spheres and characterized
by its volume, spatial boundaries, and associated residues.
The tool outputs structured pocket information suitable for downstream applications
such as docking, binding site analysis, and enzyme characterizing.
'''

# example usage:
'''
Example command:

enzywizard-pocket -i examples/input/cleaned_3GP6.cif -o examples/output/ 

'''

# input parameters:
'''
-i, --input_path
Required.
Path to the input cleaned protein structure file in CIF or PDB format.

-o, --output_dir
Required.
Directory to save the JSON report.

--min_rad
Optional.
Minimum probe radius used in PyVOL cavity detection.

Default:
  1.8

This parameter controls the smallest probe sphere used to explore cavities.
Smaller values allow detection of narrow and fine-grained pockets, but excessively
small values may lead to PyVOL failure.

--max_rad
Optional.
Maximum probe radius used in PyVOL cavity detection.

Default:
  6.2

This parameter controls the largest probe sphere used during cavity expansion.
Larger values allow identification of broader and more exposed pockets, but excessively
large values may lead to PyVOL failure.

--min_volume
Optional.
Minimum pocket volume threshold.

Default:
  50

Only pockets with volume greater than or equal to this value will be retained.
'''

# output content:
'''
The program outputs the following file into the output directory:

1. A JSON report
   - pocket_report_{protein_name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_pocket"

   - "pocket_region_statistics"
     A dictionary summarizing overall pocket properties.

     It includes:
     - "pocket_num"
       Total number of detected pockets.

     - "max_pocket_volume"
       Maximum pocket volume among all detected pockets.

     - "total_pocket_volume"
       Sum of volumes of all detected pockets.

   - "pocket_regions"
     A list describing individual pocket regions.

     Each entry contains:
     - "volume"
       Pocket volume computed by PyVOL.

     - "n_spheres"
       Number of spheres representing the pocket.

     - "residues"
       A list of residues associated with the pocket.

       Each residue contains:
       - "aa_id"
         Residue index in the cleaned structure.

       - "aa_name"
         Residue one-letter amino acid code.

     - "pocket_center_coord"
       Center coordinates of the pocket bounding box.

     - "pocket_box_boundaries"
       Size of the pocket bounding box (length, width, height).
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

3. Extract structural information
   - Retrieve the single protein chain.
   - Extract residue list and establish residue identity mapping.

4. Prepare PyVOL input
   - Convert the structure into temporary PDB file.
   - Generate a temporary PyVOL configuration file specifying probe radii and volume threshold.

5. Run PyVOL for pocket detection
   - Execute PyVOL using a rolling probe sphere algorithm.
   - Probe spheres with radii ranging from min_rad to max_rad explore the protein surface.
   - Cavities are identified based on geometric accessibility and spatial continuity.
   - Pockets are represented as clusters of overlapping spheres.

6. Parse PyVOL outputs
   - Read sphere coordinates and radii from .xyzrg files in temporary directory.
   - Extract pocket volumes from .rept report files in temporary directory.
   - Identify valid pocket object pairs (.obj and .xyzrg).

7. Compute pocket features
   - Calculate pocket center coordinates using bounding box of spheres.
   - Compute pocket bounding box dimensions.
   - Map pocket spheres to nearest residues based on spatial proximity.

8. Filter pockets
   - Remove pockets with missing volume or invalid geometry.
   - Remove pockets lacking valid residue associations.

9. Sort pockets
   - Sort valid pockets by volume in descending order.

10. Compute summary statistics
   - Calculate total number of pockets.
   - Determine maximum and total pocket volumes.

11. Save outputs
   - Generate and save a JSON report containing pocket regions and summary statistics.
'''

# dependencies:
'''
- Biopython
- PyVOL
- NumPy
'''

# references:
'''
- Guerra JV, et al. PyVOL: a PyMOL plugin for visualization, comparison,
  and volume calculation of drug-binding sites. Bioinformatics. 2021.

- PyVOL:
  https://github.com/schlessinger-lab/pyvol

- Biopython:
  https://biopython.org/
'''