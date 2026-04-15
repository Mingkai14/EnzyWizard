from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple
import json
import re

from ..utils.logging_utils import Logger
from ..utils.IO_utils import write_json_from_dict_inline_leaf_lists

from ..utils.integrate_utils import (
    SUPPORTED_OUTPUT_TYPES,
    load_json_file,
    list_json_files,
    split_integrated_graph_entries,
    validate_clean_report,
    validate_aaprops_report,
    validate_hydrocluster_report,
    validate_energy_report,
    validate_flexibility_report,
    validate_disorder_report,
    validate_conservation_report,
    validate_embedding_report,
    validate_pocket_report,
    validate_substrate_report,
    validate_dock_report,
    validate_interaction_report,
)


MUT_INTEGRATE_SUPPORTED_OUTPUT_TYPES = set(SUPPORTED_OUTPUT_TYPES) | {"enzywizard_mut_clean"}

MUT_INTEGRATE_SIDE_OUTPUT_TYPES = {
    "enzywizard_aaprops",
    "enzywizard_hydrocluster",
    "enzywizard_energy",
    "enzywizard_flexibility",
    "enzywizard_disorder",
    "enzywizard_conservation",
    "enzywizard_embedding",
    "enzywizard_pocket",
    "enzywizard_substrate",
    "enzywizard_dock",
    "enzywizard_interaction",
}


def save_mut_integrate_json(report: Dict[str, Any] | List[Any], output_path: str | Path, logger: Logger) -> bool:
    try:
        write_json_from_dict_inline_leaf_lists(report, output_path)
        return True
    except Exception as e:
        logger.print(f"[ERROR] Failed to save mut-integrate JSON: {e}")
        return False


def extract_wt_mut_protein_names_from_mutclean_report_path(
    mutclean_report_path: str | Path,
    logger: Logger,
) -> Tuple[str, str] | None:

    try:
        name = Path(mutclean_report_path).name
        m = re.fullmatch(r"mut_clean_report_(.+)\.json", name)
        if m is None:
            logger.print(
                "[ERROR] mut_clean_report file name must match "
                "mut_clean_report_{wt_protein_name}_to_{mut_protein_name}.json"
            )
            return None

        body = m.group(1).strip()
        if body == "":
            logger.print(f"[ERROR] Invalid mut_clean_report file name: {name}")
            return None

        split_token = "_to_"
        if split_token not in body:
            logger.print(
                f"[ERROR] Invalid mut_clean_report file name. Cannot find '_to_' separator: {name}"
            )
            return None

        wt_protein_name, mut_protein_name = body.split(split_token, 1)
        wt_protein_name = wt_protein_name.strip()
        mut_protein_name = mut_protein_name.strip()

        if wt_protein_name == "" or mut_protein_name == "":
            logger.print(
                f"[ERROR] Invalid WT or MUT protein name in mut_clean_report file name: {name}"
            )
            return None

        return wt_protein_name, mut_protein_name
    except Exception as e:
        logger.print(f"[ERROR] Failed to parse mut_clean_report file name: {e}")
        return None


def get_mut_integrate_supported_output_type(data: Dict[str, Any], logger: Logger) -> str | None:
    output_type = data.get("output_type")
    if not isinstance(output_type, str):
        logger.print("[ERROR] Missing or invalid output_type.")
        return None
    if output_type not in MUT_INTEGRATE_SUPPORTED_OUTPUT_TYPES:
        logger.print(f"[ERROR] Unsupported output_type for mut_integrate: {output_type}")
        return None
    return output_type


def validate_mut_integrate_report_by_type(data: Dict[str, Any], logger: Logger) -> bool:
    output_type = get_mut_integrate_supported_output_type(data, logger)
    if output_type is None:
        return False

    if output_type == "enzywizard_mut_clean":
        return validate_mutclean_report(data, logger)
    if output_type == "enzywizard_clean":
        return validate_clean_report(data, logger)
    if output_type == "enzywizard_aaprops":
        return validate_aaprops_report(data, logger)
    if output_type == "enzywizard_hydrocluster":
        return validate_hydrocluster_report(data, logger)
    if output_type == "enzywizard_energy":
        return validate_energy_report(data, logger)
    if output_type == "enzywizard_flexibility":
        return validate_flexibility_report(data, logger)
    if output_type == "enzywizard_disorder":
        return validate_disorder_report(data, logger)
    if output_type == "enzywizard_conservation":
        return validate_conservation_report(data, logger)
    if output_type == "enzywizard_embedding":
        return validate_embedding_report(data, logger)
    if output_type == "enzywizard_pocket":
        return validate_pocket_report(data, logger)
    if output_type == "enzywizard_substrate":
        return validate_substrate_report(data, logger)
    if output_type == "enzywizard_dock":
        return validate_dock_report(data, logger)
    if output_type == "enzywizard_interaction":
        return validate_interaction_report(data, logger)

    logger.print("[ERROR] Unsupported report type.")
    return False


def validate_mutclean_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_mut_clean":
        logger.print("[ERROR] mut_clean report output_type mismatch.")
        return False

    amino_acid_substitution = data.get("amino_acid_substitution")
    cleaned_amino_acid_substitution = data.get("cleaned_amino_acid_substitution")
    wt_mapping = data.get("wt_amino_acid_mapping_old_to_new")
    wt_stats = data.get("wt_clean_statistics")
    mut_mapping = data.get("mut_amino_acid_mapping_old_to_new")
    mut_stats = data.get("mut_clean_statistics")

    if not isinstance(amino_acid_substitution, str) or amino_acid_substitution.strip() == "":
        logger.print("[ERROR] Invalid amino_acid_substitution in mut_clean report.")
        return False

    if not isinstance(cleaned_amino_acid_substitution, str) or cleaned_amino_acid_substitution.strip() == "":
        logger.print("[ERROR] Invalid cleaned_amino_acid_substitution in mut_clean report.")
        return False

    if not isinstance(wt_mapping, list):
        logger.print("[ERROR] Invalid wt_amino_acid_mapping_old_to_new.")
        return False

    if not isinstance(mut_mapping, list):
        logger.print("[ERROR] Invalid mut_amino_acid_mapping_old_to_new.")
        return False

    if not isinstance(wt_stats, dict):
        logger.print("[ERROR] Invalid wt_clean_statistics.")
        return False

    if not isinstance(mut_stats, dict):
        logger.print("[ERROR] Invalid mut_clean_statistics.")
        return False

    required_stat_keys = [
        "changed_resname",
        "removed_nonstd",
        "removed_missing_bb",
        "removed_missing_heavy_atoms",
        "removed_unexpected_heavy_atoms",
        "removed_bad_occ",
        "removed_inscodes",
        "kept_residues",
    ]

    for stats, stats_name in [(wt_stats, "wt_clean_statistics"), (mut_stats, "mut_clean_statistics")]:
        for k in required_stat_keys:
            if not isinstance(stats.get(k), int):
                logger.print(f"[ERROR] Invalid {stats_name} field: {k}")
                return False

    for mapping, mapping_name in [
        (wt_mapping, "wt_amino_acid_mapping_old_to_new"),
        (mut_mapping, "mut_amino_acid_mapping_old_to_new"),
    ]:
        for item in mapping:
            if not isinstance(item, dict):
                logger.print(f"[ERROR] Invalid residue mapping entry in {mapping_name}.")
                return False

            old_residue = item.get("old_residue")
            new_residue = item.get("new_residue")

            if not isinstance(old_residue, dict) or not isinstance(new_residue, dict):
                logger.print(f"[ERROR] Invalid residue mapping structure in {mapping_name}.")
                return False

            for residue_item in [old_residue, new_residue]:
                if not isinstance(residue_item.get("aa_id"), int):
                    logger.print(f"[ERROR] Invalid aa_id in {mapping_name}.")
                    return False
                if not isinstance(residue_item.get("aa_name"), str) or residue_item["aa_name"].strip() == "":
                    logger.print(f"[ERROR] Invalid aa_name in {mapping_name}.")
                    return False
                if not isinstance(residue_item.get("hydrogen_atom_count"), int):
                    logger.print(f"[ERROR] Invalid hydrogen_atom_count in {mapping_name}.")
                    return False

    return True


def synthesize_clean_report_from_mutclean(
    mutclean_report: Dict[str, Any],
    side: str,
    logger: Logger,
) -> Dict[str, Any] | None:
    if not isinstance(mutclean_report, dict):
        logger.print("[ERROR] mut_clean_report must be a dict.")
        return None

    if side == "wt":
        mapping = mutclean_report.get("wt_amino_acid_mapping_old_to_new")
        stats = mutclean_report.get("wt_clean_statistics")
    elif side == "mut":
        mapping = mutclean_report.get("mut_amino_acid_mapping_old_to_new")
        stats = mutclean_report.get("mut_clean_statistics")
    else:
        logger.print("[ERROR] side must be 'wt' or 'mut'.")
        return None

    if not isinstance(mapping, list):
        logger.print(f"[ERROR] Missing mapping for side={side} in mut_clean_report.")
        return None

    if not isinstance(stats, dict):
        logger.print(f"[ERROR] Missing clean statistics for side={side} in mut_clean_report.")
        return None

    clean_report = {
        "output_type": "enzywizard_clean",
        "amino_acid_mapping_old_to_new": mapping,
        "clean_statistics": stats,
    }

    if not validate_clean_report(clean_report, logger):
        logger.print(f"[ERROR] Failed to synthesize clean report from mut_clean for side={side}.")
        return None

    return clean_report