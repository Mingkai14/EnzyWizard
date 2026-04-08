from __future__ import annotations
from argparse import Namespace
from ..services.dock_service import run_dock_service


def add_dock_parser(subparsers) -> None:
    parser = subparsers.add_parser("dock",help="Perform substrate docking for an input CIF/PDB file.")
    parser.add_argument("-i", "--input_path",required=True,help="Path to input CIF/PDB file.")
    parser.add_argument("-s", "--substrate_names",required=True,help='Substrate names separated by ";", for example: "lig1;lig2".')
    parser.add_argument("--substrate_dir",
        required=True,
        help="Path to a directory containing substrate SDF files."
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="Path to a directory for outputting docked structures and a JSON report."
    )
    parser.add_argument(
        "--max_docking_result_num",
        type=int,
        default=3,
        help="Maximum number of successful docking results to keep (default: 3)."
    )
    parser.add_argument(
        "--max_docking_attempt_num",
        type=int,
        default=20,
        help="Maximum number of docking attempts, where each Vina docking call counts as one attempt (default: 20)."
    )
    parser.add_argument(
        "--max_pose_read_num",
        type=int,
        default=1,
        help="Maximum number of poses to read from each successful docking run (default: 1)."
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=16,
        help="Exhaustiveness of AutoDock Vina search (default: 16)."
    )
    parser.add_argument(
        "--cpu",
        type=int,
        default=0,
        help="Number of CPUs used by AutoDock Vina (default: 0)."
    )
    parser.add_argument(
        "--min_rad",
        type=float,
        default=1.8,
        help="Minimum probe radius used in pocket detection (default: 1.8)."
    )
    parser.add_argument(
        "--max_rad",
        type=float,
        default=6.2,
        help="Maximum probe radius used in pocket detection (default: 6.2)."
    )
    parser.add_argument(
        "--min_volume",
        type=int,
        default=50,
        help="Minimum pocket volume threshold used for filtering pocket regions (default: 50)."
    )
    parser.set_defaults(func=run_dock)


def run_dock(args: Namespace) -> None:
    run_dock_service(
        input_path=args.input_path,
        substrate_names=args.substrate_names,
        substrate_dir=args.substrate_dir,
        output_dir=args.output_dir,
        max_docking_result_num=args.max_docking_result_num,
        max_docking_attempt_num=args.max_docking_attempt_num,
        max_pose_read_num=args.max_pose_read_num,
        exhaustiveness=args.exhaustiveness,
        cpu=args.cpu,
        min_rad=args.min_rad,
        max_rad=args.max_rad,
        min_volume=args.min_volume,
    )