from __future__ import annotations
import argparse

from .commands.clean import add_clean_parser
from .commands.mutclean import add_mutclean_parser
from .commands.aaprops import add_aaprops_parser
from .commands.hydrocluster import add_hydrocluster_parser
from .commands.energy import add_energy_parser
from .commands.flexibility import add_flexibility_parser
from .commands.disorder import add_disorder_parser
from .commands.conservation import add_conservation_parser
from .commands.embedding import add_embedding_parser
from .commands.pocket import add_pocket_parser
from .commands.substrate import add_substrate_parser
from .commands.dock import add_dock_parser
from .commands.interaction import add_interaction_parser

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enzywizard",
        description="EnzyWizard: an integrated toolkit for enzyme analysis."
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    add_clean_parser(subparsers)
    add_mutclean_parser(subparsers)
    add_aaprops_parser(subparsers)
    add_hydrocluster_parser(subparsers)
    add_energy_parser(subparsers)
    add_flexibility_parser(subparsers)
    add_disorder_parser(subparsers)
    add_conservation_parser(subparsers)
    add_embedding_parser(subparsers)
    add_pocket_parser(subparsers)
    add_substrate_parser(subparsers)
    add_dock_parser(subparsers)
    add_interaction_parser(subparsers)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)