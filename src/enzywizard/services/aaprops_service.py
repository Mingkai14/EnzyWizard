from __future__ import annotations
from pathlib import Path
from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists,get_stem,check_filename_length, load_protein_structure, load_dssp, write_json_from_dict_inline_leaf_lists
from ..algorithms.clean_algorithms import check_cleaned_structure
from ..algorithms.aaprops_algorithms import calculate_aa_props, calculate_aa_props_statistics, generate_aaprops_report
from ..utils.common_utils import get_optimized_filename

def run_aaprops_service(input_path: str | Path, output_dir: str | Path) ->bool:
    # ---- logger ----
    logger = Logger(output_dir)
    logger.print(f"[INFO] Aaprops processing started: {input_path}")

    # ---- check input ----
    input_path = Path(input_path)
    output_dir = Path(output_dir)

    if not file_exists(input_path):
        logger.print(f"[ERROR] Input not found: {input_path}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- get name ----
    name = get_stem(input_path)
    if not check_filename_length(name,logger):
        return False
    logger.print(f"[INFO] Protein name resolved: {name}")

    # ---- load structure ----
    structure = load_protein_structure(input_path,name,logger)
    if structure is None:
        logger.print(f"[ERROR] Failed to load structure: {input_path}")
        return False

    logger.print("[INFO] Structure loaded")

    #---- check structure ----
    if not check_cleaned_structure(structure, logger):
        return False
    logger.print(f"[INFO] Structure checked")

    #---- load dssp ----
    dssp=load_dssp(structure,logger)
    if dssp is None:
        return False

    # ---- run algorithm ----
    logger.print("[INFO] Calculating amino acid properties started")
    aa_props=calculate_aa_props(structure,dssp,logger)
    if aa_props is None:
        return False
    aa_props_statistics=calculate_aa_props_statistics(aa_props,logger)
    if aa_props_statistics is None:
        return False

    json_report_path = output_dir / get_optimized_filename(f"aaprops_report_{name}.json")

    report=generate_aaprops_report(aa_props, aa_props_statistics)
    write_json_from_dict_inline_leaf_lists(report,json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Aaprops processing finished")

    return True

