from __future__ import annotations

from argparse import Namespace
from ..services.interaction_service import run_interaction_service


def add_interaction_parser(subparsers) -> None:
    parser = subparsers.add_parser("interaction",help="Calculate protein-substrate and intra-protein interactions for an input CIF/PDB structure using docked substrate SDF files.")
    parser.add_argument("-i","--input_path",required=True,help="Path to input CIF/PDB file.")
    parser.add_argument("-s","--substrate_names",required=True,help="Input substrate names separated by ','. Each substrate name must match the corresponding docked SDF file name in substrate_dir.")
    parser.add_argument("--substrate_dir",required=True,help="Path to a directory containing docked substrate SDF files.")
    parser.add_argument("-o","--output_dir",required=True,help="Path to a directory for outputting the interaction JSON report.")
    parser.add_argument("--hbond_da_max_distance",type=float,default=3.9,help="Maximum donor-acceptor distance cutoff for hydrogen bond detection (default: 3.9).")
    parser.add_argument("--hbond_ha_max_distance",type=float,default=2.5,help="Maximum hydrogen-acceptor distance cutoff for hydrogen bond detection (default: 2.5).")
    parser.add_argument("--hbond_angle",type=float,default=90.0,help="Minimum donor-hydrogen-acceptor angle cutoff for hydrogen bond detection (default: 90.0).")
    parser.add_argument("--ionic_distance_cutoff",type=float,default=4.0,help="Maximum distance cutoff for ionic bond detection (default: 4.0).")
    parser.add_argument("--ppstack_center_distance_cutoff",type=float,default=6.5,help="Maximum ring-center distance cutoff for pi-pi stacking detection (default: 6.5).")
    parser.add_argument("--pication_distance_cutoff",type=float,default=5.0,help="Maximum ring-cation distance cutoff for pi-cation interaction detection (default: 5.0).")
    parser.add_argument("--pication_angle_cutoff",type=float,default=45.0,help="Maximum angle cutoff for pi-cation interaction detection (default: 45.0).")
    parser.add_argument("--ssbond_max_distance",type=float,default=2.5,help="Maximum sulfur-sulfur distance cutoff for disulfide bond detection (default: 2.5).")

    parser.set_defaults(func=run_interaction)


def run_interaction(args: Namespace) -> None:
    run_interaction_service(
        input_path=args.input_path,
        substrate_names=args.substrate_names,
        substrate_dir=args.substrate_dir,
        output_dir=args.output_dir,
        da_max_distance_A=args.hbond_da_max_distance,
        ha_max_distance_A=args.hbond_ha_max_distance,
        dha_min_angle_deg=args.hbond_angle,
        ionic_distance_cutoff_A=args.ionic_distance_cutoff,
        ring_center_distance_cutoff_A=args.ppstack_center_distance_cutoff,
        ring_cation_distance_cutoff_A=args.pication_distance_cutoff,
        ring_cation_angle_cutoff_deg=args.pication_angle_cutoff,
        ss_max_distance_A=args.ssbond_max_distance,
    )