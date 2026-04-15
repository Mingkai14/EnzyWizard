from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import shutil

from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists, get_stem, check_filename_length
from ..utils.common_utils import get_optimized_filename
from ..utils.mut_batch_utils import validate_mut_batch_parameter_ranges, save_mut_batch_integrate_outputs
from ..algorithms.mut_batch_algorithms import run_mut_batch_workflow


def run_mut_batch_service(
    wt_cleaned_input_path: str | Path,
    mut_cleaned_input_path: str | Path,
    wt_input_msa: str | Path,
    mut_input_msa: str | Path,
    amino_acid_substitution: str,
    substrate_names: str | None,
    wt_output_dir: str | Path,
    mut_output_dir: str | Path,
    save_extra_outputs: bool = False,
    cutoff_area: float = 10.0,
    minimize_energy: bool = True,
    minimization_iteration: int = 1000,
    energy_force_field_file: str = "charmm36.xml",
    flexibility_cutoff: float = 15.0,
    n_modes: int = 20,
    flexibility_method: str = "ANM",
    window_size: int = 11,
    min_region_length: int = 5,
    embedding_model_name: str = "esm2_t6_8M_UR50D",
    pocket_min_rad: float = 1.8,
    pocket_max_rad: float = 6.2,
    pocket_min_volume: int = 50,
    max_synonyms: int = 20,
    fp_radius: int = 2,
    n_bits: int = 512,
    num_confs: int = 5,
    prune_rms: float = 0.5,
    max_docking_attempt_num: int = 20,
    early_stop: bool = False,
    exhaustiveness: int = 16,
    cpu: int = 0,
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
    wt_cleaned_input_path = Path(wt_cleaned_input_path)
    mut_cleaned_input_path = Path(mut_cleaned_input_path)
    wt_input_msa = Path(wt_input_msa)
    mut_input_msa = Path(mut_input_msa)
    wt_output_dir = Path(wt_output_dir)
    mut_output_dir = Path(mut_output_dir)

    wt_output_dir.mkdir(parents=True, exist_ok=True)
    mut_output_dir.mkdir(parents=True, exist_ok=True)

    tmp_wt_dir_ctx: TemporaryDirectory | None = None
    tmp_mut_dir_ctx: TemporaryDirectory | None = None

    try:
        if not save_extra_outputs:
            tmp_wt_dir_ctx = TemporaryDirectory()
            tmp_mut_dir_ctx = TemporaryDirectory()
            wt_working_output_dir = Path(tmp_wt_dir_ctx.name)
            mut_working_output_dir = Path(tmp_mut_dir_ctx.name)
        else:
            wt_working_output_dir = wt_output_dir
            mut_working_output_dir = mut_output_dir

        logger = Logger(wt_working_output_dir)

        has_substrate = isinstance(substrate_names, str) and substrate_names.strip() != ""

        logger.print(
            f"[INFO] Mut_batch processing started: "
            f"wt_cleaned_input_path={wt_cleaned_input_path}, "
            f"mut_cleaned_input_path={mut_cleaned_input_path}, "
            f"wt_input_msa={wt_input_msa}, "
            f"mut_input_msa={mut_input_msa}, "
            f"substrate_names={substrate_names}, "
            f"amino_acid_substitution={amino_acid_substitution}"
        )

        if has_substrate:
            logger.print("[INFO] Substrate input detected. Full mut_batch workflow will be executed.")


        if not file_exists(wt_cleaned_input_path):
            logger.print(f"[ERROR] WT cleaned input file not found: {wt_cleaned_input_path}")
            return False

        if not file_exists(mut_cleaned_input_path):
            logger.print(f"[ERROR] MUT cleaned input file not found: {mut_cleaned_input_path}")
            return False

        if not file_exists(wt_input_msa):
            logger.print(f"[ERROR] WT input MSA file not found: {wt_input_msa}")
            return False

        if not file_exists(mut_input_msa):
            logger.print(f"[ERROR] MUT input MSA file not found: {mut_input_msa}")
            return False

        if wt_output_dir.resolve() == mut_output_dir.resolve():
            logger.print("[ERROR] wt_output_dir and mut_output_dir must be different directories.")
            return False

        if not amino_acid_substitution or not str(amino_acid_substitution).strip():
            logger.print("[ERROR] amino_acid_substitution is empty.")
            return False

        wt_protein_name = get_stem(wt_cleaned_input_path)
        if not check_filename_length(wt_protein_name, logger):
            return False

        mut_protein_name = get_stem(mut_cleaned_input_path)
        if not check_filename_length(mut_protein_name, logger):
            return False

        if wt_protein_name == mut_protein_name:
            logger.print(
                f"[ERROR] Wild-type and mutant protein names are the same: {wt_protein_name}"
            )
            return False

        wt_msa_name = get_stem(wt_input_msa)
        if not check_filename_length(wt_msa_name, logger):
            return False

        mut_msa_name = get_stem(mut_input_msa)
        if not check_filename_length(mut_msa_name, logger):
            return False

        logger.print(
            f"[INFO] Protein names resolved: wt={wt_protein_name}, mut={mut_protein_name}"
        )
        logger.print(
            f"[INFO] MSA names resolved: wt={wt_msa_name}, mut={mut_msa_name}"
        )


        if not validate_mut_batch_parameter_ranges(
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
        ):
            return False

        batch_result = run_mut_batch_workflow(
            wt_cleaned_input_path=wt_cleaned_input_path,
            mut_cleaned_input_path=mut_cleaned_input_path,
            wt_input_msa=wt_input_msa,
            mut_input_msa=mut_input_msa,
            substrate_names=substrate_names,
            amino_acid_substitution=amino_acid_substitution,
            wt_protein_name=wt_protein_name,
            mut_protein_name=mut_protein_name,
            wt_msa_name=wt_msa_name,
            mut_msa_name=mut_msa_name,
            wt_output_dir=wt_working_output_dir,
            mut_output_dir=mut_working_output_dir,
            logger=logger,
            save_extra_outputs=save_extra_outputs,
            cutoff_area=cutoff_area,
            minimize_energy=minimize_energy,
            minimization_iteration=minimization_iteration,
            energy_force_field_file=energy_force_field_file,
            flexibility_cutoff=flexibility_cutoff,
            n_modes=n_modes,
            flexibility_method=flexibility_method,
            window_size=window_size,
            min_region_length=min_region_length,
            embedding_model_name=embedding_model_name,
            pocket_min_rad=pocket_min_rad,
            pocket_max_rad=pocket_max_rad,
            pocket_min_volume=pocket_min_volume,
            max_synonyms=max_synonyms,
            fp_radius=fp_radius,
            n_bits=n_bits,
            num_confs=num_confs,
            prune_rms=prune_rms,
            max_docking_attempt_num=max_docking_attempt_num,
            early_stop=early_stop,
            exhaustiveness=exhaustiveness,
            cpu=cpu,
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
        if batch_result is None:
            return False

        mut_integrate_report = batch_result.get("mut_integrate_report")
        if not isinstance(mut_integrate_report, dict):
            logger.print("[ERROR] Missing mut_integrate_report in mut_batch result.")
            return False

        if not save_mut_batch_integrate_outputs(
                mut_integrate_report=mut_integrate_report,
                wt_output_dir=wt_working_output_dir,
                mut_output_dir=mut_working_output_dir,
                wt_protein_name=wt_protein_name,
                mut_protein_name=mut_protein_name,
                logger=logger,
        ):
            return False

        logger.print("[INFO] Mut_batch processing finished")

        if not save_extra_outputs:
            report_filename = get_optimized_filename(
                f"mut_integrate_report_{wt_protein_name}_to_{mut_protein_name}.json"
            )
            wt_nodes_filename = get_optimized_filename(f"wt_integrate_nodes_{wt_protein_name}.json")
            wt_edges_filename = get_optimized_filename(f"wt_integrate_edges_{wt_protein_name}.json")
            mut_nodes_filename = get_optimized_filename(f"mut_integrate_nodes_{mut_protein_name}.json")
            mut_edges_filename = get_optimized_filename(f"mut_integrate_edges_{mut_protein_name}.json")

            for filename in [report_filename, wt_nodes_filename, wt_edges_filename, "log.txt"]:
                src_path = wt_working_output_dir / filename
                if src_path.exists():
                    shutil.copy2(src_path, wt_output_dir / filename)

            for filename in [report_filename, mut_nodes_filename, mut_edges_filename]:
                src_path = mut_working_output_dir / filename
                if src_path.exists():
                    shutil.copy2(src_path, mut_output_dir / filename)

        wt_log_path = wt_output_dir / "log.txt"
        mut_log_path = mut_output_dir / "log.txt"

        if wt_log_path.exists():
            try:
                shutil.copy2(wt_log_path, mut_log_path)
            except Exception as e:
                logger.print(f"[ERROR] Failed to copy log.txt to mut_output_dir: {e}")
                return False

        return True

    finally:
        if tmp_wt_dir_ctx is not None:
            tmp_wt_dir_ctx.cleanup()
        if tmp_mut_dir_ctx is not None:
            tmp_mut_dir_ctx.cleanup()