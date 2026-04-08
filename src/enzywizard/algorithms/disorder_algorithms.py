from __future__ import annotations

from typing import List, Dict, Any, Tuple

from Bio.PDB.Structure import Structure

from ..utils.logging_utils import Logger
from ..utils.structure_utils import get_single_chain, get_chain_length, get_residues_by_chain, get_sequence
from ..resources.aa_physicochemical_props import hydrophobicity_dict, net_charge_dict
from ..utils.disorder_utils import moving_average

def predict_disorder_scores( sequence: str, logger: Logger, window_size: int = 11, ) -> List[float] | None:
    """
    A FoldIndex-like score:
        score = 2.785 * <hydropathy> - |<charge>| - 1.151

    score < 0  => predicted disordered
    score >= 0 => predicted ordered
    """
    hydropathy_values = [hydrophobicity_dict[aa] for aa in sequence]
    charge_values = [net_charge_dict[aa] for aa in sequence]

    mean_hydropathy = moving_average(hydropathy_values, window_size, logger)
    if mean_hydropathy is None:
        return None

    mean_charge = moving_average(charge_values, window_size, logger)
    if mean_charge is None:
        return None

    scores: List[float] = []
    for h, c in zip(mean_hydropathy, mean_charge):
        score = 2.785 * h - abs(c) - 1.151
        scores.append(score)

    return scores

def build_disordered_regions(residues: List[Tuple[Tuple[str, int, str], str, Tuple[float, float, float]]], scores: List[float], min_region_length: int, logger: Logger) -> List[Dict[str, Any]] | None:
    if len(residues) == 0:
        logger.print("[ERROR] Empty residues input in build_disordered_regions")
        return None

    if len(scores) == 0:
        logger.print("[ERROR] Empty scores input in build_disordered_regions")
        return None

    if min_region_length <= 0:
        logger.print("[ERROR] min_region_length must be positive")
        return None

    if not len(residues) == len(scores):
        logger.print("[ERROR] Length mismatch in build_disordered_regions")
        return None

    regions: List[Dict[str, Any]] = []
    start: int | None = None

    for i, score in enumerate(scores):
        is_disordered = score < 0.0

        if is_disordered and start is None:
            start = i

        elif (not is_disordered) and start is not None:
            end = i - 1

            if end - start + 1 >= min_region_length:
                residues_list: List[Dict[str, Any]] = []

                for j in range(start, end + 1):
                    key = residues[j][0]
                    resname = residues[j][1]

                    residues_list.append({
                        "aa_id": key[1],
                        "aa_name": resname,
                    })

                regions.append({
                    "length": end - start + 1,
                    "residues": residues_list,
                })

            start = None

    if start is not None:
        end = len(scores) - 1

        if end - start + 1 >= min_region_length:
            residues_list: List[Dict[str, Any]] = []

            for j in range(start, end + 1):
                key = residues[j][0]
                resname = residues[j][1]

                residues_list.append({
                    "aa_id": key[1],
                    "aa_name": resname,
                })

            regions.append({
                "length": end - start + 1,
                "residues": residues_list,
            })
    regions.sort(key=lambda x: x["length"], reverse=True)
    return regions

def compute_disordered_regions(struct: Structure,logger: Logger,window_size: int = 11,min_region_length: int = 5) -> List[Dict[str, Any]] | None:

    chain = get_single_chain(struct, logger)
    if chain is None:
        return None

    chain_length = get_chain_length(chain, logger)
    if chain_length is None:
        return None

    residues = get_residues_by_chain(chain,logger)
    if residues is None:
        return None

    if len(residues) != chain_length:
        logger.print(f"[ERROR] Sequence length mismatch in disorder calculation: chain_length={chain_length}, seq_length={len(residues)}")
        return None

    sequence = get_sequence(residues,logger)
    if sequence is None:
        return None

    scores = predict_disorder_scores(sequence=sequence,logger=logger,window_size=window_size)
    if scores is None:
        return None

    if len(scores) != len(residues):
        logger.print(f"[ERROR] Score length mismatch in disorder calculation: score_length={len(scores)}, residue_length={len(residues)}")
        return None

    regions = build_disordered_regions(
        residues=residues,
        scores=scores,
        min_region_length=min_region_length,
        logger=logger
    )

    return regions

def generate_disorder_report(disorder_regions: List[Dict[str, Any]]) -> Dict[str, Any] | None:
    return {
        "output_type": "enzywizard_disorder",
        "disorder_regions": disorder_regions,
    }