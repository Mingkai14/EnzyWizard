from __future__ import annotations
from argparse import Namespace
from ..services.aaprops_service import run_aaprops_service

def add_aaprops_parser(subparsers) -> None:
    parser = subparsers.add_parser("aaprops",help="Calculate amino acid properties from input CIF/PDB file.")
    parser.add_argument("-i","--input_path", required=True, help="Path to input CIF/PDB file.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to a directory for outputting a JSON report.")
    parser.set_defaults(func=run_aaprops)

def run_aaprops(args: Namespace) -> None:
    run_aaprops_service(input_path=args.input_path, output_dir=args.output_dir)


# input parameters:
'''
-i --input_path required input cleaned CIF/PDB protein structure file;
-o --output_dir required output directory to save JSON report
'''

# output content:
'''
The program outputs a JSON report that records:

1. "output_type": enzywizard_aaprops,

2. "aa_props": a list of amino acid–level properties for each residue,
   including:
   - amino acid index (aa_id)
   - amino acid name (one-letter code)
   - one-hot encoding of amino acid type
   - amino acid classification (multi-label and one-hot)
   - secondary structure (DSSP 8-state and one-hot)
   - relative solvent accessibility (RSA)
   - backbone dihedral angles (phi, psi)
   - physicochemical properties such as net charge, pKa, volume,
     hydrophobicity, molecular weight, and isoelectric point (pI)
   - residue coordinate (CA position)

3. "aa_props_statistics": statistical summaries over all residues,
   including:
   - counts of each amino acid type (20 categories)
   - counts of each amino acid class (8 categories, multi-label)
   - counts of each secondary structure state (DSSP 8-state)
'''

# functionality/process
'''
It processes a protein structure by:

1. Extracting a single chain from the structure;

2. Iterating over all residues in the chain and retrieving:
   - structural information from DSSP (secondary structure, RSA, phi, psi)
   - amino acid identity and classification
   - physicochemical properties from predefined dictionaries (AAindex);

3. Encoding categorical features into one-hot or multi-hot vectors;

4. Collecting all residue-level features into a list (aa_props);

5. Aggregating statistics across residues, including amino acid distribution,
   class distribution, and secondary structure distribution;

6. Generating a structured JSON report containing both detailed residue-level
   features and global statistics.
'''

# dependency:
'''
Biopython
DSSP
Predefined amino acid physicochemical properties from AAindex database
'''

# reference:
'''
Kabsch & Sander, DSSP: Dictionary of Secondary Structure of Proteins
https://swift.cmbi.umcn.nl/gv/dssp/
Biopython DSSP module documentation
https://biopython.org/docs/dev/api/Bio.PDB.DSSP.html
https://www.genome.jp/aaindex/
'''
