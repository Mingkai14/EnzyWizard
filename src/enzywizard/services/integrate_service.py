from __future__ import annotations

from pathlib import Path
from typing import Dict, Any

from ..utils.logging_utils import Logger
from ..utils.integrate_utils import (
    extract_protein_name_from_clean_report_path,
    find_unique_clean_report_path,
    get_supported_output_type,
    list_json_files,
    load_json_file,
    save_integrate_json,
    split_integrated_graph_entries,
    validate_report_by_type,
)
from ..algorithms.integrate_algorithms import integrate_reports
from ..utils.common_utils import get_optimized_filename


def run_integrate_service(input_dir: str | Path,output_dir: str | Path,strict: bool = False) -> bool:
    logger = Logger(output_dir)
    logger.print(f"[INFO] Integrate processing started: input_dir={input_dir}")

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)


    if not input_dir.exists() or not input_dir.is_dir():
        logger.print(f"[ERROR] Invalid input_dir: {input_dir}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)


    json_path_list = list_json_files(input_dir, logger)
    if json_path_list is None:
        return False

    clean_report_path = find_unique_clean_report_path(json_path_list, logger)
    if clean_report_path is None:
        return False

    protein_name = extract_protein_name_from_clean_report_path(clean_report_path, logger)
    if protein_name is None:
        return False

    unique_json_path_set = set([p.resolve() for p in json_path_list])
    if len(unique_json_path_set) > 12:
        logger.print(f"[ERROR] Number of JSON files in input_dir exceeds 12 (maximum allowed report types): {len(unique_json_path_set)}")
        return False

    report_dict: Dict[str, Dict[str, Any]] = {}

    clean_report_data = load_json_file(clean_report_path, logger)
    if clean_report_data is None:
        return False

    if not validate_report_by_type(clean_report_data, logger):
        logger.print(f"[ERROR] Invalid clean_report format: {clean_report_path}")
        return False

    clean_output_type = get_supported_output_type(clean_report_data, logger)
    if clean_output_type != "enzywizard_clean":
        logger.print(f"[ERROR] clean report in input_dir must be an enzywizard_clean JSON file: {clean_report_path}")
        return False

    report_dict["enzywizard_clean"] = clean_report_data

    for json_path in json_path_list:
        try:
            if json_path.resolve() == clean_report_path.resolve():
                continue
        except Exception:
            pass

        data = load_json_file(json_path, logger)
        if data is None:
            return False

        if not validate_report_by_type(data, logger):
            logger.print(f"[ERROR] Invalid report format: {json_path}")
            return False

        output_type = get_supported_output_type(data, logger)
        if output_type is None:
            return False

        if output_type in report_dict:
            logger.print(f"[ERROR] Duplicate report type found: {output_type}")
            return False

        report_dict[output_type] = data
        logger.print(f"[INFO] Loaded report: {json_path.name} ({output_type})")

    if strict and len(report_dict) != 12:
        logger.print(f"[ERROR] Strict mode requires exactly 12 report types, but got {len(report_dict)}.")
        return False

    integrate_report = integrate_reports(report_dict, strict, logger)
    if integrate_report is None:
        return False

    report_output_name = f"integrate_report_{protein_name}"
    nodes_output_name = f"integrate_nodes_{protein_name}"
    edges_output_name = f"integrate_edges_{protein_name}"

    json_report_path = output_dir / get_optimized_filename(f"{report_output_name}.json")
    nodes_json_path = output_dir / get_optimized_filename(f"{nodes_output_name}.json")
    edges_json_path = output_dir / get_optimized_filename(f"{edges_output_name}.json")

    if not save_integrate_json(integrate_report, json_report_path, logger):
        return False
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    integrated_graph = integrate_report.get("integrated_graph")
    if not isinstance(integrated_graph, list):
        logger.print("[ERROR] integrated_graph missing in integrate report.")
        return False

    split_result = split_integrated_graph_entries(integrated_graph, logger)
    if split_result is None:
        return False

    node_list, edge_list = split_result

    if not save_integrate_json(node_list, nodes_json_path, logger):
        return False
    logger.print(f"[INFO] Node list JSON saved: {nodes_json_path}")

    if not save_integrate_json(edge_list, edges_json_path, logger):
        return False
    logger.print(f"[INFO] Edge list JSON saved: {edges_json_path}")

    logger.print("[INFO] Integrate processing finished")
    return True