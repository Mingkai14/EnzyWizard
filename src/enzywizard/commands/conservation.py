from __future__ import annotations
from argparse import Namespace
from ..services.conservation_service import run_conservation_service

def add_conservation_parser(subparsers) -> None:
    parser = subparsers.add_parser("conservation",help="Calculate residue sequence conservation from a cleaned protein sequence, a user-provided multiple sequence alignment (MSA), and generate a detailed JSON report.")
    parser.add_argument("-i", "--input_fasta",required=True,help="Path to input cleaned protein sequence file in FASTA format.")
    parser.add_argument("-m", "--input_msa",required=True,help="Path to input multiple sequence alignment (MSA) file prepared by the user.")
    parser.add_argument("-o", "--output_dir",required=True,help="Directory to save output files, including cleaned MSA, HMM profile, and JSON report.")


    parser.set_defaults(func=run_conservation)

def run_conservation(args: Namespace) -> None:
    run_conservation_service(input_fasta=args.input_fasta, input_msa=args.input_msa, output_dir=args.output_dir)

# ==============================
# Command: enzywizard-conservation
# ==============================

# brief introduction:
'''
EnzyWizard-Conservation is a command-line tool for calculating residue
sequence conservation from a cleaned protein sequence, a user-provided
multiple sequence alignment (MSA), and generating a detailed JSON report.
It supports three formats of MSA files (Stockholm, aligned FASTA, or A3M format). 
The tool standardizes and cleans the MSA, builds a profile hidden Markov model (HMM)
using HMMER (hmmbuild), and computes conservation scores based on
HMM emission statistics. For each residue, it extracts and records: raw HMM emission log score,
normalized emission probability, and transformed Shannon entropy as a general conservation score.
'''

# example usage:
'''
Example command:

enzywizard-conservation -i examples/input/cleaned_1HVR.fasta -m examples/input/jhmmer_1HVR.sto -o examples/output/

'''

# input parameters:
'''
-i, --input_fasta
Required.
Path to input cleaned protein sequence file in FASTA format.

-m, --input_msa
Required.
Path to input multiple sequence alignment (MSA) file prepared by the user.

Supported formats:
  - Stockholm (.sto / .stockholm)
  - aligned FASTA (.fa / .fasta / .afa)
  - A3M (.a3m)

The first sequence in the MSA must match the input query sequence.

-o, --output_dir
Required.
Directory to save output files, including cleaned MSA, HMM profile, and JSON report.
'''

# output content:
'''
The program outputs the following files into the output directory:

1. A cleaned Stockholm MSA file
   - cleaned_{msa_name}.sto

2. A profile HMM file
   - hmm_profile_{msa_name}.hmm

3. A JSON report
   - conservation_report_{protein_name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_conservation"

   - "conservation_scores"
     A list describing conservation metrics for each residue.

     Each entry contains:
     - "aa_id"
       Residue index in the query sequence.

     - "aa_name"
       Residue one-letter amino acid code.

     - "hmm_emission_log_score"
       Raw emission log score from the HMM for the query amino acid.

     - "emission_probability"
       Normalized emission probability for the query amino acid.

     - "conservation_score"
       Transformed Shannon entropy computed from the full emission probability distribution.

'''

# Process:
'''
This command processes the input sequence and MSA as follows:

1. Load input sequence
   - Read the cleaned protein sequence from the FASTA file.
   - Validate sequence format and characters.

2. Load input MSA
   - Read the MSA from Stockholm, aligned FASTA, or A3M format.
   - Preserve sequence alignment structure.

3. Validate MSA
   - Confirm the MSA is non-empty and properly formatted.
   - Confirm the first sequence matches the input query sequence.
   - Validate sequence characters and alignment consistency.

4. Clean and standardize MSA
   - Normalize sequence headers.
   - Remove invalid characters.
   - Remove duplicated, empty, invalid-length, or all-gap sequences.
   - Remove lowercase insertions (for A3M).
   - Ensure all sequences are compatible with Stockholm format.

5. Save cleaned MSA
   - Write the cleaned MSA into Stockholm format.
   - Add '#=GC RF' annotation so each non-gap query residue column
     is treated as a match state in HMM construction.

6. Build HMM profile
   - Run HMMER hmmbuild with '--hand' option on the cleaned MSA.
   - Generate a profile hidden Markov model representing sequence conservation.

7. Compute conservation scores
   - Parse match emission log scores from the HMM file.
   - Convert log scores into normalized amino acid emission probabilities.
   - For each residue position:
       - extract the raw emission log score for the query amino acid
       - compute emission probability
       - compute Shannon entropy from the full probability distribution
       - convert Shannon entropy to information content

8. Save outputs
   - Generate and save a JSON report containing per-residue conservation metrics.
'''

# dependencies:
'''
- HMMER (hmmbuild)
- Biopython
- NumPy
'''

# references:
'''
- Eddy SR. Profile hidden Markov models. Bioinformatics. 1998.
- Shannon CE. A mathematical theory of communication. Bell System Technical Journal. 1948.
- HMMER:
  http://hmmer.org/
- Biopython:
  https://biopython.org/
'''