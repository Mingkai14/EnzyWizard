from __future__ import annotations

from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain

from Bio.PDB import MMCIFParser, PDBParser, Structure, Atom
from ..utils.logging_utils import Logger
from Bio.Data.IUPACData import protein_letters_3to1

from typing import List, Tuple

def get_first_model(struct: Structure, logger: Logger) -> Model | None:
    for m in struct:
        return m
    logger.print(f"[ERROR] No model found in structure")
    return None

def get_single_chain(struct: Structure, logger:Logger) -> Chain | None:
    m = get_first_model(struct, logger)
    if not m:
        return None
    for chain in m.get_chains():
        return chain
    logger.print(f"[ERROR] No chain found in structure")
    return None

def get_chain_length(chain: Chain, logger:Logger) -> int | None:
    if chain is None:
        logger.print(f"[ERROR] Bad chain input")
        return None

    length = 0

    for residue in chain:
        if residue.id[0] != " ":
            continue
        length += 1
    if not length:
        logger.print(f"[ERROR] Invalid sequence length")
        return None
    return length



def get_residues_by_chain(chain: Chain, logger: Logger) -> List[Tuple[Tuple[str,int,str], str, Tuple[float, float, float]]] | None:
    result: List[Tuple[Tuple[str, int, str], str, Tuple[float, float, float]]] = []

    for res in chain.get_residues():
        hetflag, resseq, icode = res.id

        if hetflag != " ":
            continue

        resname = res.get_resname().strip()

        # 必须有 CA 原子
        if "CA" not in res:
            logger.print(f"[ERROR] Residue {resseq} {resname} missing CA atom, skipped.")
            return None

        ca: Atom = res["CA"]
        coord = tuple(ca.get_coord())  # (x, y, z)

        result.append(((hetflag,resseq,icode), resname, coord))

    return result

def get_sequence(residues: List[Tuple[Tuple[str, int, str], str, Tuple[float, float, float]]],logger: Logger) -> str | None:

    if residues is None:
        logger.print("[ERROR] Residues input is None")
        return None

    if len(residues) == 0:
        logger.print("[ERROR] Empty residues list")
        return None

    seq_chars: List[str] = []

    for residue_key, resname, _ in residues:
        aa = protein_letters_3to1.get(resname.capitalize())
        if aa is None:
            logger.print(f"[ERROR] Unsupported residue name: {resname}")
            return None
        seq_chars.append(aa)

    return "".join(seq_chars)

