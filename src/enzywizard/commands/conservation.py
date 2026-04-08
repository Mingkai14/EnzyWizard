from __future__ import annotations
from argparse import Namespace
from ..services.conservation_service import run_conservation_service

def add_conservation_parser(subparsers) -> None:
    parser = subparsers.add_parser("conservation",help="Calculate per-residue conservation scores for an input protein sequence using its multiple sequence alignment.")
    parser.add_argument("-i", "--input_fasta",required=True,help="Path to input protein sequence file in FASTA format")
    parser.add_argument("-m", "--input_msa",required=True,help="Path to input multiple sequence alignment (MSA) file (STO/aligned FASTA/A3M format).")
    parser.add_argument("-o", "--output_dir",required=True,help="Directory to save the output JSON report, cleaned MSA file in STO format, and HMM profile file.")


    parser.set_defaults(func=run_conservation)

def run_conservation(args: Namespace) -> None:
    run_conservation_service(input_fasta=args.input_fasta, input_msa=args.input_msa, output_dir=args.output_dir)

# input parameters:
'''
-i --input_fasta Required. Path to input protein sequence file in FASTA format.

-m --input_msa Required. Path to input multiple sequence alignment (MSA) file.
Supported formats:
  - Stockholm (.sto / .stockholm)
  - aligned FASTA (.fa / .fasta / .afa)
  - A3M (.a3m)

-o --output_dir Required. Directory to save output files, including a JSON report, cleaned MSA file in STO format, and HMM profile file.
'''

# output content:
'''
The program outputs:

1. A cleaned Stockholm MSA file:
   - cleaned_{msa_name}.sto

2. A profile HMM file:
   - hmm_profile_{msa_name}.hmm

3. A JSON report:
   - conservation_report_{protein_name}.json

The JSON report contains:

- "output_type": "enzywizard_conservation"

- "conservation_scores":
  A list of per-residue conservation results.
  Each entry includes:
    - aa_id
    - aa_name
    - hmm_emission_log_score
    - emission_probability
    - conservation_score
'''

# functionality / process:
'''
This command processes a protein sequence and its MSA as follows:

1. Load input FASTA file:
   - Read the query protein sequence from the FASTA file.

2. Load input MSA file:
   - Read the MSA from Stockholm, aligned FASTA, or A3M format.

3. Validate MSA:
   - Confirm the MSA is non-empty and correctly formatted.
   - Confirm the first sequence in the MSA matches the input query sequence.
   - Confirm all sequences are valid for the corresponding MSA format.

4. Clean MSA into a unified Stockholm-compatible representation:
   - Standardize headers.
   - Remove invalid characters.
   - Remove duplicated, empty, invalid-length, or all-gap sequences.
   - Remove A3M lowercase insertion characters if present.
   - Keep only sequences compatible with a standard Stockholm alignment.

5. Save cleaned MSA:
   - Write the cleaned MSA to a Stockholm file.
   - Add '#=GC RF' annotation so that every non-gap query residue column is treated as an HMM match state.

6. Build HMM profile:
   - Run hmmbuild with '--hand' on the cleaned Stockholm file to generate a profile HMM.

7. Compute per-residue conservation scores:
   - Parse match raw emission values from the HMM file.
   - Convert emission values into normalized amino-acid probabilities.
   - For each query residue position:
       - record the HMM raw emission log score for the query amino acid
       - record the emission probability for the query amino acid
       - calculate Shannon entropy from the full emission probability distribution as the conservation score

8. Save output:
   - JSON report containing per-residue conservation results.
'''

# dependencies:
'''
- HMMER
- Biopython
'''

# reference:
'''
- Eddy SR. Profile hidden Markov models. Bioinformatics. 1998.
- Shannon CE. A mathematical theory of communication. Bell System Technical Journal. 1948.
'''