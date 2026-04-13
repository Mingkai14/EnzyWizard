from __future__ import annotations

from pathlib import Path

from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists,get_stem,check_filename_length,load_protein_structure,load_openmm_modeller,write_json_from_dict_inline_leaf_lists, load_substrate_name_and_mol_3d_list
from ..utils.common_utils import get_optimized_filename

from ..algorithms.clean_algorithms import check_cleaned_structure
from ..algorithms.interaction_algorithms import calculate_all_interaction_network,summarize_interaction_counts,generate_interaction_report


from ..utils.structure_utils import structure_has_hydrogen, structure_has_too_few_hydrogens
from ..utils.interaction_utils import filter_valid_docked_substrates


def run_interaction_service(
    input_path: str | Path,
    substrate_names: str,
    substrate_dir: str | Path,
    output_dir: str | Path,
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
    logger = Logger(output_dir)
    logger.print(f"[INFO] Interaction processing started: {input_path}")

    if bonded_h_min_distance_A < 0.5 or bonded_h_min_distance_A > 1.1:
        logger.print(f"[ERROR] bonded_h_min_distance_A out of range [0.5, 1.1]: {bonded_h_min_distance_A}")
        return False

    if bonded_h_max_distance_A < 1.0 or bonded_h_max_distance_A > 1.5:
        logger.print(f"[ERROR] bonded_h_max_distance_A out of range [1.0, 1.5]: {bonded_h_max_distance_A}")
        return False

    if bonded_h_min_distance_A >= bonded_h_max_distance_A:
        logger.print(f"[ERROR] bonded_h_min_distance_A must be < bonded_h_max_distance_A.")
        return False

    if da_max_distance_A < 2.5 or da_max_distance_A > 4.5:
        logger.print(f"[ERROR] da_max_distance_A out of range [2.5, 4.5]: {da_max_distance_A}")
        return False

    if ha_max_distance_A < 1.5 or ha_max_distance_A > 3.0:
        logger.print(f"[ERROR] ha_max_distance_A out of range [1.5, 3.0]: {ha_max_distance_A}")
        return False

    if dha_min_angle_deg < 60 or dha_min_angle_deg > 180:
        logger.print(f"[ERROR] dha_min_angle_deg out of range [60, 180]: {dha_min_angle_deg}")
        return False

    if ionic_distance_cutoff_A < 2.0 or ionic_distance_cutoff_A > 6.0:
        logger.print(f"[ERROR] ionic_distance_cutoff_A out of range [2.0, 6.0]: {ionic_distance_cutoff_A}")
        return False

    if mu < 0.001 or mu > 0.1:
        logger.print(f"[ERROR] mu out of range [0.001, 0.1]: {mu}")
        return False

    if ring_center_distance_cutoff_A < 4.0 or ring_center_distance_cutoff_A > 8.0:
        logger.print(f"[ERROR] ring_center_distance_cutoff_A out of range [4.0, 8.0]: {ring_center_distance_cutoff_A}")
        return False

    if ring_cation_distance_cutoff_A < 3.0 or ring_cation_distance_cutoff_A > 7.0:
        logger.print(f"[ERROR] ring_cation_distance_cutoff_A out of range [3.0, 7.0]: {ring_cation_distance_cutoff_A}")
        return False

    if ring_cation_angle_cutoff_deg < 0 or ring_cation_angle_cutoff_deg > 90:
        logger.print(f"[ERROR] ring_cation_angle_cutoff_deg out of range [0, 90]: {ring_cation_angle_cutoff_deg}")
        return False

    if ss_max_distance_A < 1.8 or ss_max_distance_A > 3.0:
        logger.print(f"[ERROR] ss_max_distance_A out of range [1.8, 3.0]: {ss_max_distance_A}")
        return False

    if docked_heavy_atom_distance_cutoff_A < 4.0 or docked_heavy_atom_distance_cutoff_A > 10.0:
        logger.print(
            f"[ERROR] docked_heavy_atom_distance_cutoff_A out of range [4.0, 10.0]: {docked_heavy_atom_distance_cutoff_A}")
        return False

    if min_residue_index_gap < 1 or min_residue_index_gap > 5:
        logger.print(f"[ERROR] min_residue_index_gap out of range [1, 10]: {min_residue_index_gap}")
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
        return False
    logger.print("[INFO] Structure loaded")

    modeller = load_openmm_modeller(input_path, logger)
    if modeller is None:
        return False
    logger.print("[INFO] OpenMM Modeller loaded")

    if not check_cleaned_structure(structure, logger):
        return False
    logger.print("[INFO] Structure checked")

    if not structure_has_hydrogen(structure, logger):
        logger.print("[ERROR] Protein structure does not contain hydrogen atoms. Please run 'enzywizard clean' first.")
        return False

    if structure_has_too_few_hydrogens(structure, logger):
        logger.print("[WARNING] Protein structure contains few hydrogen atoms. It is recommended to run 'enzywizard clean' first.")

    loaded = load_substrate_name_and_mol_3d_list(
        substrate_names=substrate_names,
        substrate_dir=substrate_dir,
        logger=logger,
    )
    if loaded is None:
        return False

    substrate_name_list, ligand_mol_list = loaded
    logger.print(f"[INFO] Loaded {len(ligand_mol_list)} substrate Mol(3D) object(s)")

    filtered = filter_valid_docked_substrates(
        substrate_name_list=substrate_name_list,
        ligand_mol_list=ligand_mol_list,
        modeller=modeller,
        logger=logger,
        docked_heavy_atom_distance_cutoff_A=docked_heavy_atom_distance_cutoff_A,
    )

    if filtered is None:
        return False

    valid_substrate_name_list, valid_ligand_mol_list = filtered

    logger.print(f"[INFO] Valid docked substrate count: {len(valid_ligand_mol_list)}")

    logger.print("[INFO] Interaction calculation started")
    interaction_list = calculate_all_interaction_network(
        modeller=modeller,
        ligand_mol_list=valid_ligand_mol_list,
        substrate_name_list=valid_substrate_name_list,
        struct=structure,
        logger=logger,
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
    if interaction_list is None:
        logger.print("[ERROR] Failed to calculate interaction network.")
        return False

    interaction_statistics = summarize_interaction_counts(
        interaction_list=interaction_list,
        logger=logger,
    )
    if interaction_statistics is None:
        return False


    report = generate_interaction_report(
        interaction_list=interaction_list,
        interaction_statistics=interaction_statistics,
    )

    json_name = f"interaction_report_{name}_{substrate_names}.json"
    json_name = get_optimized_filename(json_name)
    json_report_path = output_dir / json_name
    write_json_from_dict_inline_leaf_lists(report, json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Interaction processing finished")
    return True