from __future__ import annotations
from argparse import Namespace
from ..services.embedding_service import run_embedding_service

def add_embedding_parser(subparsers) -> None:
    parser = subparsers.add_parser("embedding",help="Generate protein embeddings from a cleaned protein sequence and a detailed JSON report.")
    parser.add_argument("-i", "--input_fasta",required=True,help="Path to input cleaned protein sequence file in FASTA format.")
    parser.add_argument("-o", "--output_dir",required=True,help="Directory to save the output JSON report.")
    parser.add_argument("--model_name", type=str,choices=["esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D", "esm2_t30_150M_UR50D"],default="esm2_t6_8M_UR50D",help="Model for embedding generation: esm2_t6_8M_UR50D, esm2_t12_35M_UR50D, esm2_t30_150M_UR50D.")


    parser.set_defaults(func=run_embedding)

def run_embedding(args: Namespace) -> None:
    run_embedding_service(input_fasta=args.input_fasta, output_dir=args.output_dir, model_name=args.model_name)


# ==============================
# Command: enzywizard-embedding
# ==============================

# brief introduction:
'''
EnzyWizard-Embedding is a command-line tool for generating protein
 embeddings from a cleaned protein sequence and generating a detailed JSON report.
It takes a FASTA sequence as input and uses pretrained protein language
model ESM2 to compute high-dimensional representations (embeddings) for each residue.
The tool supports three parameter scales of ESM2 model: esm2_t6_8M_UR50D, esm2_t12_35M_UR50D,
esm2_t30_150M_UR50D. The tool outputs raw floating-point values of the ESM2 embedding vectors.

'''

# example usage:
'''
Example command:

enzywizard-embedding -i examples/input/cleaned_1HVR.fasta -o examples/output/

'''

# input parameters:
'''
-i, --input_fasta
Required.
Path to input cleaned protein sequence file in FASTA format.

-o, --output_dir
Required.
Directory to save the output JSON report.

--model_name
Optional.
ESM2 model used for embedding generation.

ESM2: Evolutionary Scale Modeling 2, a transformer-based protein language model

Available options:
  - esm2_t6_8M_UR50D
  - esm2_t12_35M_UR50D
  - esm2_t30_150M_UR50D

Default:
  esm2_t6_8M_UR50D

This parameter controls the size and representational capacity of the protein language model.
Larger models provide richer representations but require more computational resources.
'''

# output content:
'''
The program outputs the following file into the output directory:

1. A JSON report
   - embedding_report_{protein_name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_embedding"

   - "embeddings"
     A list describing per-residue embedding representations.

     Embeddings: protein embeddings are numerical feature vectors that capture biochemical,
     structural, and evolutionary information learned from large-scale protein sequence data.

     Each entry contains:
     - "aa_id"
       Residue index in the input sequence.

     - "aa_name"
       Residue one-letter amino acid code.

     - "embedding"
       A list of raw floating-point values representing the ESM2 embedding vector
       for the residue.

       Embedding dimension depends on the selected model:

        - esm2_t6_8M_UR50D   (≈320-dimensional embeddings)
        - esm2_t12_35M_UR50D (≈480-dimensional embeddings)
        - esm2_t30_150M_UR50D (≈640-dimensional embeddings)


'''

# Process:
'''
This command processes the input sequence as follows:

1. Load input sequence
   - Read the cleaned protein sequence from the FASTA file.
   - Validate sequence format and ensure a single sequence is provided.

2. Validate sequence
   - Convert sequence to uppercase.
   - Ensure all residues are standard amino acids (20 canonical types).

3. Load ESM2 model
   - Load the selected pretrained ESM2 model and its alphabet.
   - ESM2 is a transformer-based protein language model trained on large-scale
     sequence databases to learn contextual residue representations.

4. Tokenize sequence
   - Convert the amino acid sequence into model-compatible tokens
     using the ESM alphabet batch converter.

5. Generate embeddings
   - Perform forward inference without gradient computation.
   - Extract hidden representations from the final transformer layer.
   - Remove special tokens and retain per-residue embeddings.

6. Format output
   - For each residue:
       - assign aa_id and aa_name
       - store the corresponding embedding vector

7. Save outputs
   - Generate and save a JSON report containing per-residue embeddings.
'''

# dependencies:
'''
- PyTorch
- fair-esm
- Biopython
'''

# references:
'''
- Lin Z, Akin H, Rao R, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science. 2023.
- Rives A, Meier J, Sercu T, et al. Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. PNAS. 2021.
- ESM (fair-esm):
  https://github.com/facebookresearch/esm
'''