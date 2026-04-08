from __future__ import annotations
from argparse import Namespace
from ..services.substrate_service import run_substrate_service


def add_substrate_parser(subparsers) -> None:
    parser = subparsers.add_parser("substrate",help="Fetch/calculate substrate molecular information and generate 3D substrate structures from substrate names/SMILES.")
    parser.add_argument("-s","--substrate_names",required=True,help="Input substrate names or SMILES strings. Multiple substrate names/SMILES should be separated by ','.")
    parser.add_argument("-o","--output_dir",required=True,help="Path to a directory for outputting a JSON report and generated substrate structure files in SDF format.")
    parser.add_argument("--max_synonyms",type=int,default=20,help="Maximum number of substrate synonyms retried when fetching SMILES from a substrate name (default: 20). A larger value may improve recall but will increase API requests and runtime.")
    parser.add_argument("--fp_radius",type=int,default=2,help="Radius used for Morgan fingerprint generation (default: 2). This controls the topological neighborhood size considered around each atom. Larger values capture broader local environments but may produce sparser fingerprints.")
    parser.add_argument("--n_bits",type=int,default=512,help="Bit size of the Morgan fingerprint vector (default: 512). Larger values reduce bit collisions but increase feature dimensionality.")
    parser.add_argument("--num_confs",type=int,default=5,help="Maximum number of 3D structures to generate for each substrate (default: 5). Larger values increase conformational coverage but also increase runtime.")
    parser.add_argument("--prune_rms",type=float,default=0.5,help="RMS threshold used to prune highly similar conformers during 3D conformer generation (default: 0.5). Smaller values keep more distinct conformers only, while larger values allow more similar conformers to be retained.")

    parser.set_defaults(func=run_substrate)


def run_substrate(args: Namespace) -> None:
    run_substrate_service(
        substrate_names=args.substrate_names,
        output_dir=args.output_dir,
        max_synonyms=args.max_synonyms,
        fp_radius=args.fp_radius,
        n_bits=args.n_bits,
        num_confs=args.num_confs,
        prune_rms=args.prune_rms,
    )


# input parameters:
'''
-s --substrate_names Required. Input substrate names or SMILES strings.
Multiple substrates are supported and should be separated by ','.

Examples:
- glucose
- CCO
- glucose,fructose,acetate
- glucose,CCO,lactate

If one input item is already a valid SMILES string, it will be recorded directly.
Its internal substrate name will be automatically assigned as smiles1, smiles2, etc.

-o --output_dir Required. Output directory for saving a JSON report and
generated substrate structure files.

--max_synonyms Optional. Maximum number of synonyms retried in matching 
when fetching a SMILES from a substrate name (default: 20) by ChEBI and PubChem APIs.
This parameter controls the upper limit of synonym-expanded retry attempts
after direct ChEBI and PubChem resolution fails. Larger values may improve
name-to-SMILES recall for difficult compound names, but they also increase
the number of web requests and runtime.

--fp_radius Optional. Radius used for Morgan fingerprint generation(default: 2).
This parameter controls how many topological bond layers around each atom are
considered when building the fingerprint. A larger radius captures broader
local structural environments.

--n_bits Optional. Bit size of the Morgan fingerprint vector (default: 512).
This parameter controls the fingerprint length. Larger values reduce bit
collisions but increase feature dimensionality.

--num_confs Optional. Maximum number of 3D conformers generated for each
substrate after hydrogen addition (default: 5).
This parameter controls conformational sampling breadth. Larger values may
capture more structural diversity but increase runtime.

--prune_rms Optional. RMS threshold used to prune highly similar conformers
during conformer embedding (default: 0.5).
This parameter controls conformer redundancy removal. Smaller values are more
strict and retain fewer, more distinct conformers; larger values allow more
similar conformers to remain.
'''


# output content:
'''
The program outputs:

1. A JSON report recording:
   - "output_type": "enzywizard_substrate"
   - "substrates": a list of substrate records

2. For each substrate record, the JSON report includes:
   - substrate_name
   - smiles
   - fingerprint
   - num_atoms
   - mol_weight
   - logp
   - structures

3. "structures" is a list of generated 3D substrate structures, where each
   entry contains:
   - structure_name
   - structure_energy

4. The actual 3D molecular objects are not stored in the JSON report.
   Instead, they are written as separate SDF files in output_dir.

5. Substrate structure file naming:
   - each generated structure is named as:
     substrate_name_1.sdf
     substrate_name_2.sdf
     substrate_name_3.sdf
     ...
   - before saving, the structure name is cleaned into a filesystem-safe name
   - invalid filename characters are replaced with '_'
   - the final structure filename stem is truncated to at most 50 characters

6. If a substrate cannot be resolved to a SMILES string, or if part of the
   downstream structure generation fails, the corresponding report fields may
   remain empty strings "" or an empty list [] depending on the step.
'''


# functionality/process
'''
It processes substrates by:

1. Parsing the input substrate_names string using ';' as the separator to
   obtain multiple substrate entries;

2. Determining whether each entry is already a valid SMILES string:
   - if yes, it is recorded directly and assigned an internal substrate name
     such as smiles1, smiles2, etc.
   - if not, it is treated as a substrate name and a SMILES string is searched

3. Resolving substrate names to SMILES strings through a staged strategy:
   - ChEBI exact search
   - exact/normalized name matching within ChEBI results
   - PubChem name-to-CID resolution
   - PubChem CID-to-SMILES extraction
   - PubChem synonym expansion followed by retrying ChEBI exact matching

4. Converting each resolved SMILES string into an RDKit 2D molecular object;

5. Calculating a Morgan fingerprint and selected 2D molecular descriptors,
   including:
   - atom count
   - molecular weight
   - logP

6. Adding explicit hydrogen atoms to the 2D molecular object;

7. Generating multiple 3D conformers from the hydrogen-added molecule using
   RDKit conformer embedding;

8. Minimizing each generated 3D conformer using the UFF force field;

9. Calculating the UFF energy of each minimized 3D conformer;

10. Saving each valid minimized 3D conformer as an individual SDF file in the
    output directory;

11. Generating a structured JSON report summarizing the resolved substrate
    information, computed descriptors, and saved 3D structure metadata.
'''


# dependency:
'''
RDKit
requests
urllib3
'''


# reference:
'''
RDKit documentation
https://www.rdkit.org/

PubChem PUG REST API
https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest

ChEBI
https://www.ebi.ac.uk/chebi/
'''