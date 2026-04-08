from __future__ import annotations
from argparse import Namespace
from ..services.embedding_service import run_embedding_service

def add_embedding_parser(subparsers) -> None:
    parser = subparsers.add_parser("embedding",help="Calculate per-residue embeddings for an input protein sequence.")
    parser.add_argument("-i", "--input_fasta",required=True,help="Path to input protein sequence file in FASTA format")
    parser.add_argument("-o", "--output_dir",required=True,help="Directory to save the output JSON report.")
    parser.add_argument("--model_name", type=str,choices=["esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D", "esm2_t30_150M_UR50D"],default="esm2_t6_8M_UR50D",help="Model for embedding generation: esm2_t6_8M_UR50D, esm2_t12_35M_UR50D, esm2_t30_150M_UR50D.")


    parser.set_defaults(func=run_embedding)

def run_embedding(args: Namespace) -> None:
    run_embedding_service(input_fasta=args.input_fasta, output_dir=args.output_dir, model_name=args.model_name)

# input parameters:
'''
-i --input_fasta Required. Path to input protein sequence file in FASTA format.

-o --output_dir Required. Directory to save the JSON report.

--model_name Optional. ESM2 model used for embedding generation.
Default: esm2_t6_8M_UR50D.
Available options:
  - esm2_t6_8M_UR50D 
  - esm2_t12_35M_UR50D 
  - esm2_t30_150M_UR50D 
This parameter controls the size and representational capacity of the protein language model.
Larger models produce longer and more informative embeddings, but require more computational resources.
'''

# output content:
'''
The program outputs:

1. A JSON report:
   - embedding_report_{name}.json

The JSON report contains:

- "output_type": "enzywizard_embedding"

- "embeddings":
  A list of per-residue embeddings.
  Each entry includes:
    - aa_id
    - aa_name
    - embedding:
        A list of floating-point values representing the learned feature vector
        of the residue (embedding dimension depends on the selected ESM2 model).
'''

# functionality / process:
'''
This command processes a protein sequence as follows:

1. Load sequence:
   - Read protein sequence from FASTA file.
   - Ensure the input contains exactly one sequence.

2. Validate sequence:
   - Convert sequence to uppercase.
   - Ensure all residues are standard amino acids (20 canonical types).

3. Load ESM2 model:
   - Load the selected pretrained ESM2 model and its alphabet.
   - Move model to the specified device.
   - Set model to evaluation mode.

4. Tokenize sequence:
   - Convert the amino acid sequence into model-compatible tokens
     using the ESM alphabet batch converter.

5. Generate embeddings:
   - Perform forward inference without gradient computation.
   - Extract representations from the final transformer layer.
   - Remove special tokens and retain per-residue embeddings.

6. Format output:
   - Map each residue to:
       - aa_id
       - aa_name 
       - embedding vector (list of floats)

7. Save output:
   - Write results into a JSON report file.
'''

# dependencies:
'''
- PyTorch
- fair-esm 
- Biopython 
'''

# reference:
'''
- Lin Z, Akin H, Rao R, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science. 2023.
'''