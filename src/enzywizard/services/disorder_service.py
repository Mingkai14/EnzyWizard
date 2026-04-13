from __future__ import annotations

from pathlib import Path

from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists,get_stem,check_filename_length,load_protein_structure, write_json_from_dict_inline_leaf_lists
from ..algorithms.clean_algorithms import check_cleaned_structure
from ..algorithms.disorder_algorithms import compute_disordered_regions, generate_disorder_report
from ..utils.common_utils import get_optimized_filename

def run_disorder_service(input_path: str | Path,output_dir: str | Path, window_size: int = 11,min_region_length: int = 5) -> bool:
    # ---- logger ----
    logger = Logger(output_dir)
    logger.print(f"[INFO] Disorder processing started: {input_path}")

    # ---- check input ----
    if not isinstance(window_size, int) or window_size < 3 or window_size > 50 or window_size % 2 == 0:
        logger.print(f"[ERROR] Invalid window_size: {window_size}. Must be odd integer in [3, 50].")
        return False

    if min_region_length < 3 or min_region_length > 50 or min_region_length > window_size:
        logger.print(f"[ERROR] Invalid min_region_length: {min_region_length}. Must be integer in [3, 50] and ≤ window_size ({window_size}).")
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
    logger.print("[INFO] Intrinsically disordered regions calculation started")
    disorder_regions=compute_disordered_regions(structure,logger,window_size=window_size,min_region_length=min_region_length)
    if disorder_regions is None:
        return False
    report = generate_disorder_report(disorder_regions, logger)
    if report is None:
        return False

    # ---- write output ----
    json_report_path = output_dir / get_optimized_filename(f"disorder_report_{name}.json")
    write_json_from_dict_inline_leaf_lists(report, json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Disorder processing finished")

    return True