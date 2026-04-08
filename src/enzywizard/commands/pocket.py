from __future__ import annotations
from argparse import Namespace
from ..services.pocket_service import run_pocket_service


def add_pocket_parser(subparsers) -> None:
    parser = subparsers.add_parser("pocket",help="Calculate pockets regions from input CIF/PDB file using PyVOL.")
    parser.add_argument("-i", "--input_path", required=True,help="Path to input CIF/PDB file.")
    parser.add_argument("-o", "--output_dir",required=True,help="Path to a directory for outputting a JSON report.")
    parser.add_argument("--min_rad",type=int,default=1.8,help="Minimum probe radius used by PyVOL for cavity detection (default: 1.8). Smaller values allow detection of narrower cavities, but excessively small values may lead to PyVOL failure.")
    parser.add_argument("--max_rad",type=int,default=6.2,help="Maximum probe radius used by PyVOL for cavity detection (default: 6.2). Larger values allow detection of broader pockets, but excessively large values may lead to PyVOL failure.")
    parser.add_argument("--min_volume",type=int,default=50,help="Minimum pocket volume threshold (default: 50). Pockets with volume below this value will be discarded.")

    parser.set_defaults(func=run_pocket)


def run_pocket(args: Namespace) -> None:
    run_pocket_service(input_path=args.input_path,output_dir=args.output_dir,min_rad=args.min_rad,max_rad=args.max_rad,min_volume=args.min_volume)

# input parameters:
'''
-i --input_path Required. Path to input protein structure file (CIF or PDB).

-o --output_dir Required. Directory to save the JSON report.

--min_rad Optional. Minimum probe radius used in PyVOL cavity detection.
Default: 1.8.
This parameter controls the smallest probe sphere used to explore cavities.
Smaller values allow detection of narrow and fine-grained pockets, but excessively small values may lead to PyVOL failure.

--max_rad Optional. Maximum probe radius used in PyVOL cavity detection.
Default: 6.2.
This parameter controls the largest probe sphere used during cavity expansion.
Larger values allow identification of broader and more exposed pockets, but excessively large values may lead to PyVOL failure.

--min_volume Optional. Minimum pocket volume threshold.
Default: 50.
Only pockets with volume greater than or equal to this value will be retained.
'''

# output content:
'''
The program outputs:

1. A JSON report:
   - pocket_report_{name}.json

The JSON report contains:

- "output_type": "enzywizard_pocket"

- "pockets":
  A list of detected pockets.
  Each entry includes:
    - volume
    - n_spheres
    - pocket_center_coord
    - pocket_box_boundaries
    - residues:
        A list of residues associated with the pocket, where each residue contains:
          - aa_id
          - aa_name
'''

# functionality / process:
'''
This command processes a protein structure as follows:

1. Load structure using Biopython.

2. Validate cleaned structure:
   - Confirm the structure contains exactly one model.
   - Confirm the structure contains exactly one chain.
   - Confirm the chain ID is "A".
   - Confirm residue numbering and formatting are valid.

3. Prepare structure for PyVOL:
   - Convert structure into PDB format.
   - Generate a temporary PyVOL configuration file.

4. Run PyVOL for pocket detection:
   - PyVOL uses a rolling probe sphere algorithm to explore cavities.
   - The protein surface is sampled using probe spheres with radii ranging
     from min_rad to max_rad.
   - Cavities are identified based on geometric accessibility and spatial continuity.
   - Pockets are represented as clusters of spheres.

5. Parse PyVOL outputs:
   - Read pocket sphere coordinates and radii from .xyzrg files.
   - Extract pocket volumes from .rept report files.
   - Count number of spheres for each pocket.

6. Compute pocket features:
   - Calculate pocket center coordinates from sphere bounding box.
   - Compute pocket bounding box size (length, width, height).
   - Map pocket spheres to nearest residues using CA atoms.

7. Filter pockets:
   - Remove pockets with missing volume or invalid geometry.
   - Remove pockets with insufficient residue mapping.

8. Sort pockets:
   - Sort all valid pockets by volume in descending order.

9. Save output:
   - JSON report containing all detected pockets and their properties.
'''

# dependencies:
'''
- Biopython
- PyVOL
- NumPy
'''

# reference:
'''
- Guerra JV, et al. PyVOL: a PyMOL plugin for visualization, comparison, and
  volume calculation of drug-binding sites. Bioinformatics. 2021.
'''