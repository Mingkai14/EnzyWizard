from __future__ import annotations
from argparse import Namespace
from ..services.aaprops_service import run_aaprops_service

def add_aaprops_parser(subparsers) -> None:
    parser = subparsers.add_parser("aaprops",help="Calculate amino acid properties from input CIF/PDB file and generate a detailed JSON report.")
    parser.add_argument("-i","--input_path", required=True, help="Path to the input cleaned protein structure file in CIF or PDB format.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to the output directory for saving the JSON report.")
    parser.set_defaults(func=run_aaprops)

def run_aaprops(args: Namespace) -> None:
    run_aaprops_service(input_path=args.input_path, output_dir=args.output_dir)


# ==============================
# Command: enzywizard-aaprops
# ==============================

# brief introduction:
'''
EnzyWizard-Aaprops is a command-line tool for calculating amino acid
properties from a cleaned protein structure and generating a detailed JSON report.
It extracts residue structural and physicochemical features, such as
amino acid identity, amino acid classification, secondary structure, relative solvent
accessibility, backbone dihedral angles, physicochemical properties predefined 
in AAindex database. It also computes overall statistics summarizing amino acid distribution,
classification distribution, and secondary structure distribution across the protein.

'''

# example usage:
'''
Example command:

enzywizard-aaprops -i examples/input/cleaned_3GP6.cif -o examples/output/

'''

# input parameters:
'''
-i, --input_path
Required.
Path to the input cleaned protein structure file in CIF or PDB format.

-o, --output_dir
Required.
Path to the output directory for saving the JSON report.
'''

# output content:
'''
The program outputs the following file into the output directory:

1. A JSON report
   - aaprops_report_{name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_aaprops"

   - "aa_props"
     A list describing amino acid properties for each residue in the
     cleaned protein structure.

     Each entry contains:
     - "aa_id"
       Residue index in the cleaned structure.

     - "aa_name"
       Residue one-letter amino acid code.

     - "aa_name_one_hot"
       One-hot encoding of the amino acid identity over the 20 standard amino acids.

     - "aa_class"
       Amino acid class label string. A residue may belong to multiple classes,
       and multiple class labels are joined by "/".

     - "aa_class_one_hot"
       Multi-hot encoding of the amino acid class membership over the predefined
       amino acid classes.

     - "aa_ss"
       DSSP 8-state secondary structure label for the residue.

     - "aa_ss_one_hot"
       One-hot encoding of the DSSP secondary structure state.

     - "aa_rsa"
       Relative solvent accessibility (RSA) value from DSSP.

     - "aa_phi"
       Backbone phi dihedral angle from DSSP.

     - "aa_psi"
       Backbone psi dihedral angle from DSSP.

     - "aa_net_charge"
       Predefined net charge value for the amino acid type.

     - "aa_pka"
       Predefined pKa value for the amino acid type.

     - "aa_volume"
       Predefined residue volume value for the amino acid type.

     - "aa_hydrophobicity"
       Predefined hydrophobicity value for the amino acid type.

     - "aa_molecular_weight"
       Predefined molecular weight value for the amino acid type.

     - "aa_pi"
       Predefined isoelectric point (pI) value for the amino acid type.

     - "aa_coord"
       Residue coordinate stored for the cleaned structure.

   - "aa_props_statistics"
     A dictionary summarizing amino acid-level statistics over the full protein.

     It includes:
     - "aa_name_statistics"
       Counts of each amino acid type over the 20 standard amino acids.

     - "aa_class_statistics"
       Counts of each predefined amino acid class across all residues.
       Because amino acid classes are multi-label, one residue may contribute
       to more than one class count.

     - "aa_ss_statistics"
       Counts of each DSSP 8-state secondary structure label.
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

3. Load DSSP
   - Run DSSP on the cleaned structure.
   - Prepare residue-level DSSP fields, including secondary structure,
     relative solvent accessibility, and backbone dihedral angles.

4. Calculate amino acid-level properties
   - Extract the single chain from the cleaned structure.
   - Iterate through all residues in chain order.
   - Convert each residue name into a normalized one-letter amino acid code.
   - Encode amino acid identity into a one-hot vector.
   - Assign amino acid class labels and encode them into a multi-hot vector.
   - Retrieve DSSP secondary structure and encode it into a one-hot vector.
   - Retrieve DSSP RSA, phi, and psi values.
   - Attach predefined physicochemical property values, including net charge,
     pKa, volume, hydrophobicity, molecular weight, and isoelectric point.
   - Store the residue coordinate together with the above features.

5. Compute summary statistics
   - Count amino acid occurrences across the 20 standard amino acid types.
   - Count amino acid class occurrences across the predefined class categories.
   - Count DSSP secondary structure occurrences across the 8-state labels.

6. Save outputs
   - Generate and save a JSON report containing both residue-level amino acid
     properties and overall summary statistics.
'''

# dependencies:
'''
- Biopython
- DSSP
- NumPy
'''

# references:
'''
- Biopython:
  https://biopython.org/

- DSSP:
  https://swift.cmbi.umcn.nl/gv/dssp/

- Biopython DSSP module documentation:
  https://biopython.org/docs/dev/api/Bio.PDB.DSSP.html

- AAindex:
  https://www.genome.jp/aaindex/
'''
