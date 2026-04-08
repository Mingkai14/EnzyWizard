from __future__ import annotations
from argparse import Namespace
from ..services.disorder_service import run_disorder_service


def add_disorder_parser(subparsers) -> None:
    parser = subparsers.add_parser("disorder",help="Calculate intrinsically disordered regions from input CIF/PDB file.")
    parser.add_argument("-i", "--input_path",required=True,help="Path to input CIF/PDB file.")
    parser.add_argument("-o", "--output_dir",required=True,help="Path to a directory for outputting a JSON report.")
    parser.add_argument("--window_size",type=int,default=11,help="Sliding window size for FoldIndex-like disorder score calculation (default: 11). Larger values produce smoother regional trends, while smaller values are more sensitive to local fluctuations.")
    parser.add_argument("--min_region_length",type=int,default=5,help="Minimum number of consecutive residues required to define a disordered region (default: 5). Shorter predicted segments will be ignored.")

    parser.set_defaults(func=run_disorder)


def run_disorder(args: Namespace) -> None:
    run_disorder_service(input_path=args.input_path,output_dir=args.output_dir,window_size=args.window_size,min_region_length=args.min_region_length)



# input parameters:
'''
-i --input_path Required. Path to input protein structure file (CIF or PDB).

-o --output_dir Required. Directory to save the JSON report.

--window_size Optional. Sliding window size used for FoldIndex-like score smoothing.
Default: 11.
This parameter controls how many neighboring residues are used when averaging
hydrophobicity and net charge. Larger values emphasize broad regional trends,
while smaller values make the prediction more sensitive to local variation.

--min_region_length Optional. Minimum length of a predicted disordered segment.
Default: 5.
Only consecutive residue segments with at least this many residues and
FoldIndex-like scores below zero are reported as disordered regions.
'''

# output content:
'''
The program outputs:

1. A JSON report:
   - disorder_report_{name}.json

The JSON report contains:

- "output_type": "enzywizard_disorder"

- "disordered_regions":
  A list of predicted intrinsically disordered regions.
  Each entry includes:
    - length
    - residues:
        A list of residues in the region, where each residue contains:
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
   - Confirm residue numbering and formatting are valid for downstream analysis.

3. Extract residue sequence:
   - Collect standard protein residues from the structure.
   - Convert residue names from three-letter codes to one-letter sequence.

4. Predict intrinsic disorder using a FoldIndex-like algorithm:
   - Assign each residue a hydrophobicity value and a net charge value (from AAindex database).
   - Compute sliding-window average hydrophobicity and net charge.
   - Calculate a FoldIndex-like score for each residue:
       score = 2.785 * <mean hydrophobicity> - |<mean net charge>| - 1.151
   - Interpret residues with score < 0 as disorder-prone.

5. Build disordered regions:
   - Merge consecutive disorder-prone residues into continuous regions.
   - Keep only regions whose length is at least min_region_length.

6. Save output:
   - JSON report containing all predicted disordered regions.
'''

# dependencies:
'''
- Biopython
- AAindex
'''

# reference:
'''
- Prilusky J, Felder CE, Zeev-Ben-Mordehai T, et al. FoldIndex©: a simple tool to predict whether a given protein sequence is intrinsically unfolded. Bioinformatics. 2005.
'''