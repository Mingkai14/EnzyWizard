from __future__ import annotations

from pathlib import Path

from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists,get_stem,check_filename_length,load_protein_structure,write_json_from_dict_inline_leaf_lists

from ..algorithms.clean_algorithms import check_cleaned_structure
from ..algorithms.dock_algorithms import dock_multiple_substrates_from_structure,save_docking_results_and_generate_dock_report
from ..utils.common_utils import get_optimized_filename


def run_dock_service(
    input_path: str | Path,
    substrate_names: str,
    substrate_dir: str | Path,
    output_dir: str | Path,
    max_docking_attempt_num: int = 20,
    early_stop: bool = False,
    exhaustiveness: int = 16,
    cpu: int = 0,
    min_rad: float = 1.8,
    max_rad: float = 6.2,
    min_volume: int = 50,
) -> bool:
    logger = Logger(output_dir)
    logger.print(f"[INFO] Dock processing started: {input_path}")

    if max_docking_attempt_num <= 0 or max_docking_attempt_num > 100 or exhaustiveness <= 0 or exhaustiveness > 64 or min_rad < 1.2 or min_volume <= 20 or min_rad >= max_rad:
        logger.print(
            f"[ERROR] Invalid docking parameters. Require: max_docking_attempt_num (1–100), exhaustiveness (1–64), min_rad ≥ 1.2, max_rad > min_rad, min_volume > 20.")
        return False

    input_path = Path(input_path)
    substrate_dir = Path(substrate_dir)
    output_dir = Path(output_dir)

    if not file_exists(input_path):
        logger.print(f"[ERROR] Input not found: {input_path}")
        return False

    if not substrate_names or not str(substrate_names).strip():
        logger.print("[ERROR] substrate_names is empty.")
        return False

    if not substrate_dir.exists() or not substrate_dir.is_dir():
        logger.print(f"[ERROR] Invalid substrate_dir: {substrate_dir}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    name = get_stem(input_path)
    if not check_filename_length(name, logger):
        return False
    logger.print(f"[INFO] Protein name resolved: {name}")

    structure = load_protein_structure(input_path, name, logger)
    if structure is None:
        logger.print(f"[ERROR] Failed to load structure: {input_path}")
        return False
    logger.print("[INFO] Structure loaded")

    if not check_cleaned_structure(structure, logger):
        return False
    logger.print("[INFO] Structure checked")

    logger.print("[INFO] Docking workflow started")
    docking_result_list = dock_multiple_substrates_from_structure(
        struct=structure,
        substrate_names=substrate_names,
        substrate_dir=substrate_dir,
        logger=logger,
        max_docking_attempt_num=max_docking_attempt_num,
        early_stop=early_stop,
        exhaustiveness=exhaustiveness,
        cpu=cpu,
        min_rad=min_rad,
        max_rad=max_rad,
        min_volume=min_volume,
    )
    if docking_result_list is None:
        return False

    logger.print(f"[INFO] Docking finished")

    logger.print("[INFO] Saving docking results and generating report")
    report = save_docking_results_and_generate_dock_report(
        docking_result_list=docking_result_list,
        struct=structure,
        protein_name=name,
        output_dir=output_dir,
        logger=logger,
    )
    if report is None:
        return False

    json_name = f"dock_report_{name}_{substrate_names}.json"
    json_name = get_optimized_filename(json_name)
    json_report_path = output_dir / json_name
    write_json_from_dict_inline_leaf_lists(report, json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Dock processing finished")
    return True