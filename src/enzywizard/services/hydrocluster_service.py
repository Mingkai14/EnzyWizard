from __future__ import annotations

from pathlib import Path

from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists,get_stem,check_filename_length,load_protein_structure,write_json_from_dict_inline_leaf_lists
from ..algorithms.clean_algorithms import check_cleaned_structure

from ..algorithms.hydrocluster_algorithms import compute_hydrophobic_clusters,generate_hydrocluster_report
from ..utils.common_utils import get_optimized_filename


def run_hydrocluster_service(input_path: str | Path,output_dir: str | Path, cutoff_area: float = 10.0) -> bool:
    # ---- logger ----
    logger = Logger(output_dir)
    logger.print(f"[INFO] Hydrocluster processing started: {input_path}")

    # ---- check input ----
    if cutoff_area <= 0:
        logger.print(f"[ERROR] Invalid cutoff_area: {cutoff_area}. Must be a positive number.")
        return False

    input_path = Path(input_path)
    output_dir = Path(output_dir)

    if not file_exists(input_path):
        logger.print(f"[ERROR] Input not found: {input_path}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- get name ----
    name = get_stem(input_path)
    if not check_filename_length(name, logger):
        return False
    logger.print(f"[INFO] Protein name resolved: {name}")

    # ---- load structure ----
    structure = load_protein_structure(input_path, name, logger)
    if structure is None:
        logger.print(f"[ERROR] Failed to load structure: {input_path}")
        return False

    logger.print("[INFO] Structure loaded")

    #---- check structure ----
    if not check_cleaned_structure(structure, logger):
        return False
    logger.print(f"[INFO] Structure checked")

    # ---- run algorithm ----
    logger.print("[INFO] Hydrophobic cluster calculation started")
    clusters = compute_hydrophobic_clusters(struct=structure,logger=logger,cutoff_area=cutoff_area)
    if clusters is None:
        return False

    # ---- generate report ----
    report = generate_hydrocluster_report(clusters=clusters,struct=structure,logger=logger)
    if report is None:
        return False

    # ---- write output ----
    json_report_path = output_dir / get_optimized_filename(f"hydrocluster_report_{name}.json")
    write_json_from_dict_inline_leaf_lists(report, json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Hydrocluster processing finished")

    return True