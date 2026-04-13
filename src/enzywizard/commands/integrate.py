from __future__ import annotations
from argparse import Namespace
from ..services.integrate_service import run_integrate_service

def add_integrate_parser(subparsers) -> None:
    parser = subparsers.add_parser("integrate",help="Integrate EnzyWizard JSON reports.")
    parser.add_argument("-i", "--clean_report_path",required=True,help="Path to the required clean report JSON file.")
    parser.add_argument("-d", "--input_dir",required=True,help="Path to a directory containing JSON reports to integrate.")
    parser.add_argument("-o", "--output_dir",required=True,help="Path to output directory for integrated JSON and GraphML.")
    parser.add_argument("--strict",type=lambda x: str(x).lower() in ["true", "1", "yes"],default=False,help="Whether to require all 12 report types and all node fields (True/False, default: False).")
    parser.set_defaults(func=run_integrate)

def run_integrate(args: Namespace) -> None:
    run_integrate_service(clean_report_path=args.clean_report_path, input_dir=args.input_dir, output_dir=args.output_dir, strict=args.strict)