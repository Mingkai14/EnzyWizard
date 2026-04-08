from __future__ import annotations

from pathlib import Path
from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists,get_stem,check_filename_length,load_protein_structure, write_json_from_dict_inline_leaf_lists
from ..algorithms.clean_algorithms import check_cleaned_structure
from ..algorithms.pocket_algorithms import compute_pockets, generate_pocket_report

def run_pocket_service(input_path: str | Path, output_dir: str | Path, min_rad: int = 1.8, max_rad: int =6.2, min_volume: int =50) -> bool:
    # ---- logger ----
    logger = Logger(output_dir)
    logger.print(f"[INFO] Pocket processing started: {input_path}")

    # ---- check input ----
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
    logger.print("[INFO] Pocket regions calculation started")
    pocket_regions=compute_pockets(structure,logger,min_rad=min_rad,max_rad=max_rad,min_volume=min_volume)
    if pocket_regions is None:
        return False
    report=generate_pocket_report(pocket_regions)

    # ---- write output ----
    json_report_path = output_dir / f"pocket_report_{name}.json"
    write_json_from_dict_inline_leaf_lists(report, json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Pocket processing finished")

    return True