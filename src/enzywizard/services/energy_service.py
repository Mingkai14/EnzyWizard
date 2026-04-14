from __future__ import annotations

from pathlib import Path
from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists,get_stem,check_filename_length,load_protein_structure,load_openmm_structure,write_json_from_dict_inline_leaf_lists
from ..algorithms.clean_algorithms import check_cleaned_structure
from ..algorithms.energy_algorithms import compute_energy_terms, generate_energy_report
from ..utils.common_utils import get_optimized_filename

def run_energy_service(input_path: str | Path,output_dir: str | Path, minimize_energy: bool = True, minimization_iteration:int = 1000,force_field_file="charmm36.xml") -> bool:
    # ---- logger ----
    logger = Logger(output_dir)
    logger.print(f"[INFO] Energy processing started: {input_path}")

    # ---- check input ----
    if minimization_iteration <= 0:
        logger.print(f"[ERROR] Invalid minimization_iteration: {minimization_iteration}. Must be a positive integer.")
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

    openmm_structure = load_openmm_structure(input_path, logger)
    if openmm_structure is None:
        logger.print(f"[ERROR] Failed to load structure: {input_path}")
        return False

    logger.print("[INFO] Structure loaded")

    #---- check structure ----
    if not check_cleaned_structure(structure, logger):
        return False
    logger.print(f"[INFO] Structure checked")

    # ---- run algorithm ----
    logger.print("[INFO] Energy calculation started")
    energy_terms = compute_energy_terms(struct=openmm_structure,logger=logger,minimize_energy=minimize_energy,minimization_iteration=minimization_iteration,force_field_file=force_field_file)
    if energy_terms is None:
        return False

    # ---- generate report ----
    report = generate_energy_report(energy_terms=energy_terms,logger=logger)
    if report is None:
        return False

    # ---- write output ----
    json_report_path = output_dir / get_optimized_filename(f"energy_report_{name}.json")
    write_json_from_dict_inline_leaf_lists(report, json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Energy processing finished")

    return True

