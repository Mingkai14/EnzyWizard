from __future__ import annotations

from pathlib import Path
from typing import Dict, Any
import shutil

from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists
from ..utils.common_utils import get_optimized_filename

from ..utils.mut_integrate_utils import (
    load_json_file,
    list_json_files,
    save_mut_integrate_json,
    split_integrated_graph_entries,
    extract_wt_mut_protein_names_from_mutclean_report_path,
    validate_mut_integrate_report_by_type,
    get_mut_integrate_supported_output_type,
    MUT_INTEGRATE_SIDE_OUTPUT_TYPES,
)

from ..algorithms.mut_integrate_algorithms import integrate_mut_reports


def run_mut_integrate_service(mutclean_report_path: str | Path,wt_input_dir: str | Path,mut_input_dir: str | Path,wt_output_dir: str | Path,mut_output_dir: str | Path,strict: bool = False) -> bool:
    wt_output_dir = Path(wt_output_dir)
    mut_output_dir = Path(mut_output_dir)

    wt_output_dir.mkdir(parents=True, exist_ok=True)
    mut_output_dir.mkdir(parents=True, exist_ok=True)

    logger = Logger(wt_output_dir)
    logger.print(
        f"[INFO] Mut_integrate processing started: "
        f"mut_clean_report={mutclean_report_path}, "
        f"wt_input_dir={wt_input_dir}, "
        f"mut_input_dir={mut_input_dir}, "
        f"wt_output_dir={wt_output_dir}, "
        f"mut_output_dir={mut_output_dir}"
    )
    mutclean_report_path = Path(mutclean_report_path)
    wt_input_dir = Path(wt_input_dir)
    mut_input_dir = Path(mut_input_dir)

    if not file_exists(mutclean_report_path):
        logger.print(f"[ERROR] mut_clean_report not found: {mutclean_report_path}")
        return False

    if not wt_input_dir.exists() or not wt_input_dir.is_dir():
        logger.print(f"[ERROR] Invalid wt_input_dir: {wt_input_dir}")
        return False

    if not mut_input_dir.exists() or not mut_input_dir.is_dir():
        logger.print(f"[ERROR] Invalid mut_input_dir: {mut_input_dir}")
        return False

    if wt_output_dir.resolve() == mut_output_dir.resolve():
        logger.print("[ERROR] wt_output_dir and mut_output_dir must be different directories.")
        return False

    protein_name_pair = extract_wt_mut_protein_names_from_mutclean_report_path(mutclean_report_path, logger)
    if protein_name_pair is None:
        return False
    wt_protein_name, mut_protein_name = protein_name_pair

    mutclean_report_data = load_json_file(mutclean_report_path, logger)
    if mutclean_report_data is None:
        return False

    if not validate_mut_integrate_report_by_type(mutclean_report_data, logger):
        logger.print(f"[ERROR] Invalid mut_clean_report format: {mutclean_report_path}")
        return False

    mutclean_output_type = get_mut_integrate_supported_output_type(mutclean_report_data, logger)
    if mutclean_output_type != "enzywizard_mut_clean":
        logger.print(
            f"[ERROR] mut_clean_report_path must point to an enzywizard_mut_clean JSON file: {mutclean_report_path}"
        )
        return False

    wt_json_path_list = list_json_files(wt_input_dir, logger)
    if wt_json_path_list is None:
        return False

    mut_json_path_list = list_json_files(mut_input_dir, logger)
    if mut_json_path_list is None:
        return False

    wt_unique_json_path_set = set([p.resolve() for p in wt_json_path_list])
    mut_unique_json_path_set = set([p.resolve() for p in mut_json_path_list])

    if len(wt_unique_json_path_set) > 12:
        logger.print(f"[ERROR] Number of JSON files in wt_input_dir exceeds 12 (maximum allowed report types): {len(wt_unique_json_path_set)}")
        return False

    if len(mut_unique_json_path_set) > 12:
        logger.print(f"[ERROR] Number of JSON files in mut_input_dir exceeds 12 (maximum allowed report types): {len(mut_unique_json_path_set)}")
        return False

    wt_report_dict: Dict[str, Dict[str, Any]] = {}
    mut_report_dict: Dict[str, Dict[str, Any]] = {}

    for json_path in wt_json_path_list:
        if json_path.resolve() == mutclean_report_path.resolve():
            continue

        data = load_json_file(json_path, logger)
        if data is None:
            return False

        if not validate_mut_integrate_report_by_type(data, logger):
            logger.print(f"[ERROR] Invalid report format in wt_input_dir: {json_path}")
            return False

        output_type = get_mut_integrate_supported_output_type(data, logger)
        if output_type is None:
            return False

        if output_type not in MUT_INTEGRATE_SIDE_OUTPUT_TYPES:
            logger.print(
                f"[ERROR] Unsupported report type in wt_input_dir for mut_integrate: {output_type} ({json_path.name})")
            return False

        if output_type in wt_report_dict:
            logger.print(f"[ERROR] Duplicate report type found in wt_input_dir: {output_type}")
            return False

        wt_report_dict[output_type] = data
        logger.print(f"[INFO] Loaded wt report: {json_path.name} ({output_type})")

    for json_path in mut_json_path_list:
        if json_path.resolve() == mutclean_report_path.resolve():
            continue

        data = load_json_file(json_path, logger)
        if data is None:
            return False

        if not validate_mut_integrate_report_by_type(data, logger):
            logger.print(f"[ERROR] Invalid report format in mut_input_dir: {json_path}")
            return False

        output_type = get_mut_integrate_supported_output_type(data, logger)
        if output_type is None:
            return False

        if output_type not in MUT_INTEGRATE_SIDE_OUTPUT_TYPES:
            logger.print(
                f"[ERROR] Unsupported report type in mut_input_dir for mut_integrate: {output_type} ({json_path.name})")
            return False

        if output_type in mut_report_dict:
            logger.print(f"[ERROR] Duplicate report type found in mut_input_dir: {output_type}")
            return False

        mut_report_dict[output_type] = data
        logger.print(f"[INFO] Loaded mut report: {json_path.name} ({output_type})")

    if strict:
        expected_side_count = len(MUT_INTEGRATE_SIDE_OUTPUT_TYPES)
        if len(wt_report_dict) != expected_side_count:
            logger.print(f"[ERROR] Strict mode requires exactly {expected_side_count} report types in wt_input_dir, but got {len(wt_report_dict)}.")
            return False
        if len(mut_report_dict) != expected_side_count:
            logger.print(f"[ERROR] Strict mode requires exactly {expected_side_count} report types in mut_input_dir, but got {len(mut_report_dict)}.")
            return False

    mut_integrate_report = integrate_mut_reports(
        mutclean_report=mutclean_report_data,
        wt_report_dict=wt_report_dict,
        mut_report_dict=mut_report_dict,
        strict=strict,
        logger=logger,
    )
    if mut_integrate_report is None:
        return False

    report_output_name = f"mut_integrate_report_{wt_protein_name}_to_{mut_protein_name}"
    wt_nodes_output_name = f"wt_integrate_nodes_{wt_protein_name}"
    wt_edges_output_name = f"wt_integrate_edges_{wt_protein_name}"
    mut_nodes_output_name = f"mut_integrate_nodes_{mut_protein_name}"
    mut_edges_output_name = f"mut_integrate_edges_{mut_protein_name}"

    wt_report_json_path = wt_output_dir / get_optimized_filename(f"{report_output_name}.json")
    mut_report_json_path = mut_output_dir / get_optimized_filename(f"{report_output_name}.json")

    wt_nodes_json_path = wt_output_dir / get_optimized_filename(f"{wt_nodes_output_name}.json")
    wt_edges_json_path = wt_output_dir / get_optimized_filename(f"{wt_edges_output_name}.json")

    mut_nodes_json_path = mut_output_dir / get_optimized_filename(f"{mut_nodes_output_name}.json")
    mut_edges_json_path = mut_output_dir / get_optimized_filename(f"{mut_edges_output_name}.json")

    if not save_mut_integrate_json(mut_integrate_report, wt_report_json_path, logger):
        return False
    logger.print(f"[INFO] Mut_integrate report JSON saved: {wt_report_json_path}")

    if not save_mut_integrate_json(mut_integrate_report, mut_report_json_path, logger):
        return False
    logger.print(f"[INFO] Mut_integrate report JSON saved: {mut_report_json_path}")

    wt_integrated_graph = mut_integrate_report.get("wt_integrated_graph")
    if not isinstance(wt_integrated_graph, list):
        logger.print("[ERROR] wt_integrated_graph missing in mut-integrate report.")
        return False

    mut_integrated_graph = mut_integrate_report.get("mut_integrated_graph")
    if not isinstance(mut_integrated_graph, list):
        logger.print("[ERROR] mut_integrated_graph missing in mut-integrate report.")
        return False

    wt_split_result = split_integrated_graph_entries(wt_integrated_graph, logger)
    if wt_split_result is None:
        return False
    wt_node_list, wt_edge_list = wt_split_result

    mut_split_result = split_integrated_graph_entries(mut_integrated_graph, logger)
    if mut_split_result is None:
        return False
    mut_node_list, mut_edge_list = mut_split_result

    if not save_mut_integrate_json(wt_node_list, wt_nodes_json_path, logger):
        return False
    logger.print(f"[INFO] WT node list JSON saved: {wt_nodes_json_path}")

    if not save_mut_integrate_json(wt_edge_list, wt_edges_json_path, logger):
        return False
    logger.print(f"[INFO] WT edge list JSON saved: {wt_edges_json_path}")

    if not save_mut_integrate_json(mut_node_list, mut_nodes_json_path, logger):
        return False
    logger.print(f"[INFO] MUT node list JSON saved: {mut_nodes_json_path}")

    if not save_mut_integrate_json(mut_edge_list, mut_edges_json_path, logger):
        return False
    logger.print(f"[INFO] MUT edge list JSON saved: {mut_edges_json_path}")

    logger.print("[INFO] Mut_integrate processing finished")

    wt_log_path = wt_output_dir / "log.txt"
    mut_log_path = mut_output_dir / "log.txt"

    if wt_log_path.exists():
        try:
            shutil.copy2(wt_log_path, mut_log_path)
        except Exception as e:
            logger.print(f"[ERROR] Failed to copy log.txt to mut_output_dir: {e}")
            return False

    return True