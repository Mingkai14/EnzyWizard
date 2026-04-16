from __future__ import annotations
from argparse import Namespace
from ..services.disorder_service import run_disorder_service


def add_disorder_parser(subparsers) -> None:
    parser = subparsers.add_parser("disorder",help="Predict intrinsically disordered regions from a cleaned protein structure and generating a detailed JSON report.")
    parser.add_argument("-i", "--input_path",required=True,help="Path to the input cleaned protein structure file in CIF or PDB format.")
    parser.add_argument("-o", "--output_dir",required=True,help="Path to the output directory for saving the JSON report.")
    parser.add_argument("--window_size",type=int,default=11,help="Sliding window size for disorder score calculation (default: 11). Larger values produce smoother regional trends, while smaller values are more sensitive to local fluctuations.")
    parser.add_argument("--min_region_length",type=int,default=5,help="Minimum number of consecutive residues required to define a disordered region (default: 5). Shorter predicted segments will be ignored.")

    parser.set_defaults(func=run_disorder)


def run_disorder(args: Namespace) -> None:
    run_disorder_service(input_path=args.input_path,output_dir=args.output_dir,window_size=args.window_size,min_region_length=args.min_region_length)



# ==============================
# Command: enzywizard-disorder
# ==============================

# brief introduction:
'''
EnzyWizard-Disorder is a command-line tool for predicting intrinsically
disordered regions from a cleaned protein structure and generating a detailed
JSON report. It extracts the protein sequence from the cleaned structure and
applies a FoldIndex-like disorder prediction algorithm based on residue
hydrophobicity and net charge. The program identifies continuous residue
segments with disorder-prone physicochemical signatures and summarizes both
region-level predictions and overall disorder statistics across the protein.

'''

# example usage:
'''
Example command:

enzywizard-disorder -i examples/input/cleaned_3GP6.cif -o examples/output/

'''

# input parameters:
'''
-i, --input_path
Required.
Path to the input cleaned protein structure file in CIF or PDB format.

-o, --output_dir
Required.
Path to the output directory for saving the JSON report.

--window_size
Optional.
Sliding window size for disorder score calculation.
Default: 11.
This parameter controls how many neighboring residues are used when averaging
hydrophobicity and net charge along the sequence. Larger values produce smoother
regional trends and emphasize broad disorder tendencies, while smaller values
are more sensitive to local sequence fluctuations.

--min_region_length
Optional.
Minimum number of consecutive residues required to define a disordered region.
Default: 5.
Only continuous residue segments whose predicted disorder scores remain below
zero and whose lengths are at least this threshold are reported as disordered
regions.
'''

# output content:
'''
The program outputs the following file into the output directory:

1. A JSON report
   - disorder_report_{name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_disorder"

   - "disorder_region_statistics"
     A dictionary summarizing disorder-region-level statistics over the full protein.

     It includes:
     - "region_num"
       Total number of predicted disordered regions.

     - "max_region_length"
       Length of the longest predicted disordered region.

     - "total_region_length"
       Total number of residues covered by all predicted disordered regions.

   - "disorder_regions"
     A list describing predicted intrinsically disordered regions in the
     cleaned protein structure.

     Each entry contains:
     - "length"
       Number of residues in the predicted disordered region.

     - "residues"
       A list of residues belonging to the region.

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
   - Validate that the input structure satisfies the cleaned-structure requirement.

3. Extract sequence information
   - Extract the single chain from the cleaned structure.
   - Retrieve residues in chain order from the cleaned protein structure.
   - Confirm consistency between chain length and residue list length.
   - Convert residue names into a normalized one-letter amino acid sequence.

4. Predict residue-level disorder propensity
   - Assign each residue a predefined hydrophobicity value and net charge value.
   - Apply a sliding-window moving average to hydrophobicity and net charge
     along the sequence.
   - Calculate a FoldIndex-like score for each residue using:

       score = 2.785 × <mean hydrophobicity> - |<mean net charge>| - 1.151

   - Interpret the score as follows:
     - score < 0: predicted disorder-prone residue
     - score >= 0: predicted ordered residue

   - This algorithm is based on the principle that intrinsically disordered
     protein regions tend to have relatively low hydrophobicity and relatively
     high net charge. Lower hydrophobicity weakens the driving force for stable
     compact folding, while higher net charge increases electrostatic repulsion,
     both of which favor structural disorder.

5. Build disordered regions
   - Scan the residue-level FoldIndex-like scores in sequence order.
   - Merge consecutive residues with scores below zero into continuous
     candidate disordered segments.
   - Remove segments shorter than the user-defined minimum region length.
   - Record the residue composition of each retained disordered region.
   - Sort predicted regions by descending region length.

6. Compute disorder statistics
   - Count the total number of predicted disordered regions.
   - Calculate the maximum predicted region length.
   - Calculate the total number of residues covered by predicted disordered regions.

7. Save outputs
   - Generate and save a JSON report containing both predicted disordered
     regions and overall disorder-region statistics.
'''

# dependencies:
'''
- Biopython
- AAindex
'''

# references:
'''
- FoldIndex:
  Prilusky J, Felder CE, Zeev-Ben-Mordehai T, et al.
  FoldIndex©: a simple tool to predict whether a given protein sequence is
  intrinsically unfolded.
  Bioinformatics. 2005.

- AAindex:
  https://www.genome.jp/aaindex/

- Biopython:
  https://biopython.org/
'''