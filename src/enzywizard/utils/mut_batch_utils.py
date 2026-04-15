from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict, Any

from ..utils.logging_utils import Logger
from ..utils.batch_utils import validate_batch_parameter_ranges
from ..utils.common_utils import get_optimized_filename
from ..utils.mut_integrate_utils import (
    save_mut_integrate_json,
    split_integrated_graph_entries,
)


def validate_mut_batch_parameter_ranges(
    logger: Logger,
    cutoff_area: float = 10.0,
    minimize_energy: bool = True,
    minimization_iteration: int = 1000,
    flexibility_cutoff: float = 15.0,
    n_modes: int = 20,
    window_size: int = 11,
    min_region_length: int = 5,
    pocket_min_rad: float = 1.8,
    pocket_max_rad: float = 6.2,
    pocket_min_volume: int = 50,
    max_synonyms: int = 20,
    fp_radius: int = 2,
    n_bits: int = 512,
    num_confs: int = 5,
    prune_rms: float = 0.5,
    max_docking_attempt_num: int = 20,
    exhaustiveness: int = 16,
    dock_min_rad: float = 1.8,
    dock_max_rad: float = 6.2,
    dock_min_volume: int = 50,
    bonded_h_min_distance_A: float = 0.8,
    bonded_h_max_distance_A: float = 1.3,
    da_max_distance_A: float = 3.9,
    ha_max_distance_A: float = 2.5,
    dha_min_angle_deg: float = 90.0,
    ionic_distance_cutoff_A: float = 4.0,
    mu: float = 0.01,
    ring_center_distance_cutoff_A: float = 6.5,
    ring_cation_distance_cutoff_A: float = 5.0,
    ring_cation_angle_cutoff_deg: float = 45.0,
    ss_max_distance_A: float = 2.5,
    docked_heavy_atom_distance_cutoff_A: float = 6.5,
    min_residue_index_gap: int = 3,
) -> bool:
    return validate_batch_parameter_ranges(
        logger=logger,
        cutoff_area=cutoff_area,
        minimize_energy=minimize_energy,
        minimization_iteration=minimization_iteration,
        flexibility_cutoff=flexibility_cutoff,
        n_modes=n_modes,
        window_size=window_size,
        min_region_length=min_region_length,
        pocket_min_rad=pocket_min_rad,
        pocket_max_rad=pocket_max_rad,
        pocket_min_volume=pocket_min_volume,
        max_synonyms=max_synonyms,
        fp_radius=fp_radius,
        n_bits=n_bits,
        num_confs=num_confs,
        prune_rms=prune_rms,
        max_docking_attempt_num=max_docking_attempt_num,
        exhaustiveness=exhaustiveness,
        dock_min_rad=dock_min_rad,
        dock_max_rad=dock_max_rad,
        dock_min_volume=dock_min_volume,
        bonded_h_min_distance_A=bonded_h_min_distance_A,
        bonded_h_max_distance_A=bonded_h_max_distance_A,
        da_max_distance_A=da_max_distance_A,
        ha_max_distance_A=ha_max_distance_A,
        dha_min_angle_deg=dha_min_angle_deg,
        ionic_distance_cutoff_A=ionic_distance_cutoff_A,
        mu=mu,
        ring_center_distance_cutoff_A=ring_center_distance_cutoff_A,
        ring_cation_distance_cutoff_A=ring_cation_distance_cutoff_A,
        ring_cation_angle_cutoff_deg=ring_cation_angle_cutoff_deg,
        ss_max_distance_A=ss_max_distance_A,
        docked_heavy_atom_distance_cutoff_A=docked_heavy_atom_distance_cutoff_A,
        min_residue_index_gap=min_residue_index_gap,
    )


def copy_substrate_sdf_files(
    source_dir: str | Path,
    target_dir: str | Path,
    resolved_substrate_names: str,
    logger: Logger,
) -> bool:
    source_dir = Path(source_dir)
    target_dir = Path(target_dir)

    if not source_dir.exists() or not source_dir.is_dir():
        logger.print(f"[ERROR] Invalid substrate source_dir: {source_dir}")
        return False

    substrate_name_list = [x.strip() for x in resolved_substrate_names.split(",") if x.strip()]
    if len(substrate_name_list) == 0:
        logger.print("[ERROR] substrate_names is empty after parsing.")
        return False

    target_dir.mkdir(parents=True, exist_ok=True)

    copied_path_set = set()
    copied_count = 0

    for substrate_name in substrate_name_list:
        if not isinstance(substrate_name, str) or substrate_name.strip() == "":
            logger.print("[ERROR] Invalid substrate_name in substrate_names.")
            return False
        substrate_name = substrate_name.strip()

        matched_files = []

        single_path = source_dir / f"{substrate_name}.sdf"
        if single_path.exists() and single_path.is_file():
            matched_files.append(single_path)

        matched_files.extend(sorted(source_dir.glob(f"{substrate_name}_*.sdf")))

        unique_matched_files = []
        seen_resolved = set()
        for path in matched_files:
            resolved = path.resolve()
            if resolved in seen_resolved:
                continue
            seen_resolved.add(resolved)
            unique_matched_files.append(path)

        if len(unique_matched_files) == 0:
            logger.print(f"[ERROR] No substrate SDF files found for substrate: {substrate_name}")
            return False

        for sdf_path in unique_matched_files:
            target_path = target_dir / sdf_path.name
            try:
                shutil.copy2(sdf_path, target_path)
            except Exception as e:
                logger.print(f"[ERROR] Failed to copy substrate SDF file {sdf_path} -> {target_dir}: {e}")
                return False

            copied_key = target_path.resolve()
            if copied_key not in copied_path_set:
                copied_path_set.add(copied_key)
                copied_count += 1

    logger.print(f"[INFO] Copied {copied_count} substrate SDF file(s) to {target_dir}")
    return True



def save_mut_batch_integrate_outputs(
    mut_integrate_report: Dict[str, Any],
    wt_output_dir: str | Path,
    mut_output_dir: str | Path,
    wt_protein_name: str,
    mut_protein_name: str,
    logger: Logger,
) -> bool:
    wt_output_dir = Path(wt_output_dir)
    mut_output_dir = Path(mut_output_dir)

    report_filename = get_optimized_filename(
        f"mut_integrate_report_{wt_protein_name}_to_{mut_protein_name}.json"
    )
    wt_nodes_filename = get_optimized_filename(f"wt_integrate_nodes_{wt_protein_name}.json")
    wt_edges_filename = get_optimized_filename(f"wt_integrate_edges_{wt_protein_name}.json")
    mut_nodes_filename = get_optimized_filename(f"mut_integrate_nodes_{mut_protein_name}.json")
    mut_edges_filename = get_optimized_filename(f"mut_integrate_edges_{mut_protein_name}.json")

    wt_report_path = wt_output_dir / report_filename
    mut_report_path = mut_output_dir / report_filename

    wt_nodes_path = wt_output_dir / wt_nodes_filename
    wt_edges_path = wt_output_dir / wt_edges_filename

    mut_nodes_path = mut_output_dir / mut_nodes_filename
    mut_edges_path = mut_output_dir / mut_edges_filename

    if not save_mut_integrate_json(mut_integrate_report, wt_report_path, logger):
        return False
    logger.print(f"[INFO] Mut_integrate report JSON saved: {wt_report_path}")

    if not save_mut_integrate_json(mut_integrate_report, mut_report_path, logger):
        return False
    logger.print(f"[INFO] Mut_integrate report JSON saved: {mut_report_path}")

    wt_integrated_graph = mut_integrate_report.get("wt_integrated_graph")
    if not isinstance(wt_integrated_graph, list):
        logger.print("[ERROR] wt_integrated_graph missing in mut_batch integrate report.")
        return False

    mut_integrated_graph = mut_integrate_report.get("mut_integrated_graph")
    if not isinstance(mut_integrated_graph, list):
        logger.print("[ERROR] mut_integrated_graph missing in mut_batch integrate report.")
        return False

    wt_split_result = split_integrated_graph_entries(wt_integrated_graph, logger)
    if wt_split_result is None:
        return False
    wt_node_list, wt_edge_list = wt_split_result

    mut_split_result = split_integrated_graph_entries(mut_integrated_graph, logger)
    if mut_split_result is None:
        return False
    mut_node_list, mut_edge_list = mut_split_result

    if not save_mut_integrate_json(wt_node_list, wt_nodes_path, logger):
        return False
    logger.print(f"[INFO] WT node list JSON saved: {wt_nodes_path}")

    if not save_mut_integrate_json(wt_edge_list, wt_edges_path, logger):
        return False
    logger.print(f"[INFO] WT edge list JSON saved: {wt_edges_path}")

    if not save_mut_integrate_json(mut_node_list, mut_nodes_path, logger):
        return False
    logger.print(f"[INFO] MUT node list JSON saved: {mut_nodes_path}")

    if not save_mut_integrate_json(mut_edge_list, mut_edges_path, logger):
        return False
    logger.print(f"[INFO] MUT edge list JSON saved: {mut_edges_path}")

    return True