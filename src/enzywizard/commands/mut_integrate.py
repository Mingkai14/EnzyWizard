from __future__ import annotations
from argparse import Namespace
from ..services.mut_integrate_service import run_mut_integrate_service


def add_mut_integrate_parser(subparsers) -> None:
    parser = subparsers.add_parser("mut_integrate",help="Integrate wild-type and mutant EnzyWizard JSON reports.",)
    parser.add_argument("-i","--mutclean_report_path",required=True,help="Path to the required mutclean report JSON file.",)
    parser.add_argument("-w","--wt_input_dir",required=True,help="Path to a directory containing wild-type JSON reports to integrate.",)
    parser.add_argument("-m","--mut_input_dir",required=True,help="Path to a directory containing mutant JSON reports to integrate.",)
    parser.add_argument("-o","--output_dir",required=True,help="Path to output directory for mut-integrated JSON files.",)
    parser.add_argument("--strict",type=lambda x: str(x).lower() in ["true", "1", "yes"],default=False,help="Whether to require all report types on both wild-type and mutant sides (True/False, default: False).",)
    parser.set_defaults(func=run_mut_integrate)


def run_mut_integrate(args: Namespace) -> None:
    run_mut_integrate_service(
        mutclean_report_path=args.mutclean_report_path,
        wt_input_dir=args.wt_input_dir,
        mut_input_dir=args.mut_input_dir,
        output_dir=args.output_dir,
        strict=args.strict,
    )