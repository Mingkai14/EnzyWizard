from __future__ import annotations

from pathlib import Path
from typing import Dict, Any

from openmm.app import Modeller
from Bio.PDB.Structure import Structure

from ..utils.logging_utils import Logger
from ..utils.IO_utils import (
    load_protein_structure,
    load_dssp,
    load_msa,
    write_msa,
    write_hmm,
    save_substrate_structures,
    load_sdf_mol_3d,
    structure_to_pdbfile,
)
from ..utils.batch_utils import build_batch_output_paths
from ..utils.mut_batch_utils import (
    copy_substrate_sdf_files
)
from ..utils.structure_utils import (
    structure_has_hydrogen,
    get_fasta_dict_from_structure,
    get_single_chain,
    get_chain_length,
)
from ..utils.sequence_utils import check_msa, clean_msa_to_sto
from ..utils.interaction_utils import filter_valid_docked_substrates
from ..utils.substrate_utils import build_docked_mol_from_atom_info
from ..utils.mut_clean_utils import check_amino_acid_substitution
from ..algorithms.mut_clean_algorithms import (
    get_cleaned_amino_acid_substitution,
    generate_mutclean_report,
)

from ..algorithms.clean_algorithms import (
    clean_structure_to_single_chain_A,
    generate_clean_report,
    check_cleaned_structure,
    validate_clean_mapping_coordinates,
)
from ..algorithms.aaprops_algorithms import (
    calculate_aa_props,
    calculate_aa_props_statistics,
    generate_aaprops_report,
)
from ..algorithms.hydrocluster_algorithms import (
    compute_hydrophobic_clusters,
    generate_hydrocluster_report,
)
from ..algorithms.energy_algorithms import compute_energy_terms, generate_energy_report
from ..algorithms.flexibility_algorithms import (
    compute_protein_rmsf,
    generate_flexibility_report,
)
from ..algorithms.disorder_algorithms import (
    compute_disordered_regions,
    generate_disorder_report,
)
from ..algorithms.conservation_algorithms import (
    compute_conservation_scores,
    generate_conservation_report,
)
from ..algorithms.embedding_algorithms import (
    generate_embedding,
    generate_embedding_report,
)
from ..algorithms.pocket_algorithms import compute_pockets, generate_pocket_report
from ..algorithms.substrate_algorithms import (
    get_substrate_dict_list_from_input,
    get_completed_smiles_list,
    get_substrate_feature_list,
    generate_substrate_report,
)
from ..algorithms.dock_algorithms import (
    dock_multiple_substrates_from_structure,
    save_docking_results_and_generate_dock_report,
)
from ..algorithms.interaction_algorithms import (
    calculate_all_interaction_network,
    summarize_interaction_counts,
    generate_interaction_report,
)
from ..algorithms.mut_integrate_algorithms import integrate_mut_reports


def _run_mut_batch_side_workflow(
    cleaned_structure: Structure,
    clean_report: Dict[str, Any],
    input_msa: str | Path,
    protein_name: str,
    msa_name: str,
    output_dir: str | Path,
    substrate_names: str | None,
    substrate_report: Dict[str, Any] | None,
    substrate_dir: str | Path | None,
    logger: Logger,
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
) -> Dict[str, Dict[str, Any]] | None:
    input_msa = Path(input_msa)
    output_dir = Path(output_dir)
    substrate_dir = Path(substrate_dir) if substrate_dir is not None else None

    report_dict: Dict[str, Dict[str, Any]] = {}
    report_dict["enzywizard_clean"] = clean_report

    has_substrate = (
            isinstance(substrate_names, str)
            and substrate_names.strip() != ""
            and isinstance(substrate_report, dict)
            and substrate_dir is not None
    )

    if not check_cleaned_structure(cleaned_structure, logger):
        logger.print(f"[ERROR] Cleaned structure failed validation: {protein_name}")
        return None

    if not structure_has_hydrogen(cleaned_structure, logger):
        logger.print(f"[ERROR] Cleaned structure does not contain hydrogen atoms: {protein_name}")
        return None

    cleaned_pdbfile = structure_to_pdbfile(cleaned_structure, logger, protein_name=protein_name)
    if cleaned_pdbfile is None:
        return None

    try:
        cleaned_modeller = Modeller(cleaned_pdbfile.topology, cleaned_pdbfile.positions)
    except Exception as e:
        logger.print(f"[ERROR] Failed to build OpenMM Modeller for {protein_name}: {e}")
        return None
    logger.print(f"[INFO] OpenMM Modeller loaded: {protein_name}")

    sequence_dict = get_fasta_dict_from_structure(cleaned_structure, logger, header=protein_name)
    if sequence_dict is None:
        return None
    logger.print(f"[INFO] Sequence prepared: {protein_name}")

    logger.print(f"[INFO] Aaprops calculation started: {protein_name}")
    dssp = load_dssp(cleaned_structure, logger)
    if dssp is None:
        return None

    aa_props = calculate_aa_props(cleaned_structure, dssp, logger)
    if aa_props is None:
        return None

    aa_props_statistics = calculate_aa_props_statistics(aa_props, logger)
    if aa_props_statistics is None:
        return None

    aaprops_report = generate_aaprops_report(aa_props, aa_props_statistics)
    report_dict["enzywizard_aaprops"] = aaprops_report

    logger.print(f"[INFO] Hydrocluster calculation started: {protein_name}")
    clusters = compute_hydrophobic_clusters(cleaned_structure, logger, cutoff_area=cutoff_area)
    if clusters is None:
        return None

    hydrocluster_report = generate_hydrocluster_report(clusters, cleaned_structure, logger)
    if hydrocluster_report is None:
        return None
    report_dict["enzywizard_hydrocluster"] = hydrocluster_report

    logger.print(f"[INFO] Energy calculation started: {protein_name}")
    cleaned_pdbfile = structure_to_pdbfile(cleaned_structure, logger, protein_name=protein_name)
    if cleaned_pdbfile is None:
        return None

    energy_terms = compute_energy_terms(
        struct=cleaned_pdbfile,
        logger=logger,
        minimize_energy=minimize_energy,
        minimization_iteration=minimization_iteration,
        force_field_file=energy_force_field_file,
    )
    if energy_terms is None:
        return None

    energy_report = generate_energy_report(energy_terms=energy_terms, logger=logger)
    if energy_report is None:
        return None
    report_dict["enzywizard_energy"] = energy_report

    logger.print(f"[INFO] Flexibility calculation started: {protein_name}")
    protein_rmsf = compute_protein_rmsf(
        cleaned_structure,
        logger,
        cutoff=flexibility_cutoff,
        n_modes=n_modes,
        method=flexibility_method,
    )
    if protein_rmsf is None:
        return None

    flexibility_report = generate_flexibility_report(protein_rmsf)
    report_dict["enzywizard_flexibility"] = flexibility_report

    logger.print(f"[INFO] Disorder calculation started: {protein_name}")
    disorder_regions = compute_disordered_regions(
        cleaned_structure,
        logger,
        window_size=window_size,
        min_region_length=min_region_length,
    )
    if disorder_regions is None:
        return None

    disorder_report = generate_disorder_report(disorder_regions, logger)
    if disorder_report is None:
        return None
    report_dict["enzywizard_disorder"] = disorder_report

    logger.print(f"[INFO] Conservation calculation started: {protein_name}")
    msa_list = load_msa(input_msa, logger)
    if msa_list is None:
        return None

    if not check_msa(input_msa, sequence_dict, msa_list, logger):
        return None

    cleaned_msa_list = clean_msa_to_sto(msa_list, logger)
    if cleaned_msa_list is None:
        return None

    path_dict = build_batch_output_paths(protein_name=protein_name, msa_name=msa_name, output_dir=output_dir)

    if not write_msa(cleaned_msa_list, path_dict["cleaned_sto"], logger):
        return None
    logger.print(f"[INFO] Cleaned MSA STO file saved: {path_dict['cleaned_sto']}")

    if not write_hmm(path_dict["cleaned_sto"], path_dict["hmm"], logger):
        return None
    logger.print(f"[INFO] HMM Profile file saved: {path_dict['hmm']}")

    conservation_scores = compute_conservation_scores(path_dict["hmm"], sequence_dict, logger)
    if conservation_scores is None:
        return None

    conservation_report = generate_conservation_report(conservation_scores)
    report_dict["enzywizard_conservation"] = conservation_report

    logger.print(f"[INFO] Embedding calculation started: {protein_name}")
    embeddings = generate_embedding(sequence_dict, logger, model_name=embedding_model_name)
    if embeddings is None:
        return None

    embedding_report = generate_embedding_report(embeddings)
    report_dict["enzywizard_embedding"] = embedding_report

    logger.print(f"[INFO] Pocket calculation started: {protein_name}")
    pocket_regions = compute_pockets(
        cleaned_structure,
        logger,
        min_rad=pocket_min_rad,
        max_rad=pocket_max_rad,
        min_volume=pocket_min_volume,
    )
    if pocket_regions is None:
        return None

    pocket_report = generate_pocket_report(pocket_regions)
    report_dict["enzywizard_pocket"] = pocket_report

    if has_substrate:
        report_dict["enzywizard_substrate"] = substrate_report


    if has_substrate:
        logger.print(f"[INFO] Docking workflow started: {protein_name}")
        docking_result_list = dock_multiple_substrates_from_structure(
            struct=cleaned_structure,
            substrate_names=substrate_names,
            substrate_dir=substrate_dir,
            logger=logger,
            max_docking_attempt_num=max_docking_attempt_num,
            early_stop=early_stop,
            exhaustiveness=exhaustiveness,
            cpu=cpu,
            min_rad=dock_min_rad,
            max_rad=dock_max_rad,
            min_volume=dock_min_volume,
        )
        if docking_result_list is None:
            return None

        dock_report = save_docking_results_and_generate_dock_report(
            docking_result_list=docking_result_list,
            struct=cleaned_structure,
            protein_name=protein_name,
            output_dir=output_dir,
            logger=logger,
        )
        if dock_report is None:
            return None
        report_dict["enzywizard_dock"] = dock_report

        if len(docking_result_list) == 0:
            logger.print(f"[ERROR] Empty docking_result_list for interaction workflow: {protein_name}")
            return None

        ligand_mol_list = []
        substrate_name_list = []

        docking_result = docking_result_list[0]

        for ligand in docking_result["docked_substrate_info_list"]:
            substrate_name = ligand["substrate_name"]
            atom_info_list = ligand["atom_info_list"]
            source_sdf_path = ligand["source_sdf_path"]

            original_mol = load_sdf_mol_3d(source_sdf_path, logger)
            if original_mol is None:
                return None

            docked_mol = build_docked_mol_from_atom_info(
                original_mol,
                atom_info_list,
                logger,
            )
            if docked_mol is None:
                return None

            ligand_mol_list.append(docked_mol)
            substrate_name_list.append(substrate_name)

        logger.print(f"[INFO] Loaded {len(ligand_mol_list)} docked substrate Mol(3D) object(s): {protein_name}")

        filtered = filter_valid_docked_substrates(
            substrate_name_list=substrate_name_list,
            ligand_mol_list=ligand_mol_list,
            modeller=cleaned_modeller,
            logger=logger,
            docked_heavy_atom_distance_cutoff_A=docked_heavy_atom_distance_cutoff_A,
        )
        if filtered is None:
            return None

        valid_substrate_name_list, valid_ligand_mol_list = filtered
        logger.print(f"[INFO] Valid docked substrate count: {len(valid_ligand_mol_list)} ({protein_name})")
    else:
        valid_substrate_name_list = []
        valid_ligand_mol_list = []
        logger.print(
            f"[INFO] No substrate input detected. Only intra-protein interactions will be calculated: {protein_name}")

    logger.print(f"[INFO] Interaction workflow started: {protein_name}")

    interaction_list = calculate_all_interaction_network(
        modeller=cleaned_modeller,
        ligand_mol_list=valid_ligand_mol_list,
        substrate_name_list=valid_substrate_name_list,
        struct=cleaned_structure,
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
        logger.print(f"[ERROR] Failed to calculate interaction network: {protein_name}")
        return None

    interaction_statistics = summarize_interaction_counts(interaction_list=interaction_list, logger=logger)
    if interaction_statistics is None:
        return None

    interaction_report = generate_interaction_report(
        interaction_list=interaction_list,
        interaction_statistics=interaction_statistics,
    )
    report_dict["enzywizard_interaction"] = interaction_report

    logger.print(f"[INFO] Side workflow finished: {protein_name}")
    return report_dict


def run_mut_batch_workflow(
    wt_cleaned_input_path: str | Path,
    mut_cleaned_input_path: str | Path,
    wt_input_msa: str | Path,
    mut_input_msa: str | Path,
    substrate_names: str | None,
    amino_acid_substitution: str,
    wt_protein_name: str,
    mut_protein_name: str,
    wt_msa_name: str,
    mut_msa_name: str,
    wt_output_dir: str | Path,
    mut_output_dir: str | Path,
    logger: Logger,
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
) -> Dict[str, Any] | None:
    wt_cleaned_input_path = Path(wt_cleaned_input_path)
    mut_cleaned_input_path = Path(mut_cleaned_input_path)
    wt_input_msa = Path(wt_input_msa)
    mut_input_msa = Path(mut_input_msa)
    wt_output_dir = Path(wt_output_dir)
    mut_output_dir = Path(mut_output_dir)
    has_substrate = isinstance(substrate_names, str) and substrate_names.strip() != ""

    wt_output_dir.mkdir(parents=True, exist_ok=True)
    mut_output_dir.mkdir(parents=True, exist_ok=True)

    logger.print("[INFO] Mut_batch workflow started")

    try:
        wt_structure = load_protein_structure(wt_cleaned_input_path, wt_protein_name, logger)
    except Exception:
        wt_structure = None

    try:
        mut_structure = load_protein_structure(mut_cleaned_input_path, mut_protein_name, logger)
    except Exception:
        mut_structure = None

    if wt_structure is None:
        logger.print(f"[ERROR] Failed to load WT structure: {wt_cleaned_input_path}")
        return None

    if mut_structure is None:
        logger.print(f"[ERROR] Failed to load MUT structure: {mut_cleaned_input_path}")
        return None

    if not check_cleaned_structure(wt_structure, logger):
        logger.print("[ERROR] WT input structure is not a valid cleaned structure.")
        return None

    if not check_cleaned_structure(mut_structure, logger):
        logger.print("[ERROR] MUT input structure is not a valid cleaned structure.")
        return None

    if not structure_has_hydrogen(wt_structure, logger):
        logger.print("[ERROR] WT input cleaned structure does not contain hydrogen atoms.")
        return None

    if not structure_has_hydrogen(mut_structure, logger):
        logger.print("[ERROR] MUT input cleaned structure does not contain hydrogen atoms.")
        return None

    wt_chain = get_single_chain(wt_structure, logger)
    mut_chain = get_single_chain(mut_structure, logger)
    if wt_chain is None or mut_chain is None:
        return None

    wt_seq_length = get_chain_length(wt_chain, logger)
    mut_seq_length = get_chain_length(mut_chain, logger)
    if wt_seq_length is None or mut_seq_length is None:
        return None

    if wt_seq_length != mut_seq_length:
        logger.print(
            f"[ERROR] WT and MUT sequence lengths are not equal: {wt_seq_length} vs {mut_seq_length}"
        )
        return None

    if not check_amino_acid_substitution(
        amino_acid_substitution,
        wt_length=wt_seq_length,
        mut_length=mut_seq_length,
        logger=logger,
    ):
        return None

    wt_clean_result = clean_structure_to_single_chain_A(wt_structure, logger)
    if wt_clean_result is None:
        return None
    wt_cleaned_structure, wt_mapping_old_to_new, wt_clean_stats = wt_clean_result

    mut_clean_result = clean_structure_to_single_chain_A(mut_structure, logger)
    if mut_clean_result is None:
        return None
    mut_cleaned_structure, mut_mapping_old_to_new, mut_clean_stats = mut_clean_result

    if not check_cleaned_structure(wt_cleaned_structure, logger):
        return None

    if not check_cleaned_structure(mut_cleaned_structure, logger):
        return None

    if not validate_clean_mapping_coordinates(
        wt_structure, wt_cleaned_structure, wt_mapping_old_to_new, logger
    ):
        return None

    if not validate_clean_mapping_coordinates(
        mut_structure, mut_cleaned_structure, mut_mapping_old_to_new, logger
    ):
        return None

    if not structure_has_hydrogen(wt_cleaned_structure, logger):
        logger.print("[ERROR] WT structure does not contain hydrogen atoms.")
        return None

    if not structure_has_hydrogen(mut_cleaned_structure, logger):
        logger.print("[ERROR] MUT structure does not contain hydrogen atoms.")
        return None

    wt_cleaned_chain = get_single_chain(wt_cleaned_structure, logger)
    mut_cleaned_chain = get_single_chain(mut_cleaned_structure, logger)
    if wt_cleaned_chain is None or mut_cleaned_chain is None:
        return None

    wt_cleaned_length = get_chain_length(wt_cleaned_chain, logger)
    mut_cleaned_length = get_chain_length(mut_cleaned_chain, logger)
    if wt_cleaned_length is None or mut_cleaned_length is None:
        return None

    if wt_cleaned_length != mut_cleaned_length:
        logger.print(
            f"[ERROR] Cleaned WT and MUT sequence lengths are not equal: "
            f"{wt_cleaned_length} vs {mut_cleaned_length}"
        )
        return None

    cleaned_amino_acid_substitution = get_cleaned_amino_acid_substitution(
        wt_mapping_old_to_new=wt_mapping_old_to_new,
        mut_mapping_old_to_new=mut_mapping_old_to_new,
        mutation=amino_acid_substitution,
        logger=logger,
    )
    if cleaned_amino_acid_substitution is None:
        return None

    wt_clean_report = generate_clean_report(
        wt_structure,
        wt_cleaned_structure,
        wt_mapping_old_to_new,
        wt_clean_stats,
        logger,
    )
    if wt_clean_report is None:
        return None

    mut_clean_report = generate_clean_report(
        mut_structure,
        mut_cleaned_structure,
        mut_mapping_old_to_new,
        mut_clean_stats,
        logger,
    )
    if mut_clean_report is None:
        return None

    mutclean_report = generate_mutclean_report(
        old_aas=amino_acid_substitution,
        cleaned_aas=cleaned_amino_acid_substitution,
        wt_structure=wt_structure,
        wt_cleaned_structure=wt_cleaned_structure,
        mut_structure=mut_structure,
        mut_cleaned_structure=mut_cleaned_structure,
        wt_mapping_old_to_new=wt_mapping_old_to_new,
        wt_stats=wt_clean_stats,
        mut_mapping_old_to_new=mut_mapping_old_to_new,
        mut_stats=mut_clean_stats,
        logger=logger,
    )
    if mutclean_report is None:
        return None

    resolved_substrate_names: str | None = None
    substrate_report: Dict[str, Any] | None = None

    if has_substrate:
        logger.print("[INFO] Substrate generation started")
        substrate_dict_list = get_substrate_dict_list_from_input(substrate_names, logger)
        if substrate_dict_list is None:
            return None

        substrate_dict_list = get_completed_smiles_list(
            substrate_dict_list,
            logger,
            max_synonyms=max_synonyms,
        )
        if substrate_dict_list is None:
            return None

        substrate_feature_list = get_substrate_feature_list(
            substrate_dict_list,
            logger,
            fp_radius=fp_radius,
            n_bits=n_bits,
            num_confs=num_confs,
            prune_rms=prune_rms,
        )
        if substrate_feature_list is None:
            return None

        resolved_substrate_names = ",".join(item["substrate_name"] for item in substrate_dict_list)

        if not save_substrate_structures(substrate_feature_list, wt_output_dir, logger):
            return None
        logger.print(f"[INFO] Substrate structures saved to WT side: {wt_output_dir}")

        if not copy_substrate_sdf_files(wt_output_dir, mut_output_dir, resolved_substrate_names, logger):
            return None
        logger.print(f"[INFO] Substrate structures copied to MUT side: {mut_output_dir}")

        substrate_report = generate_substrate_report(substrate_feature_list, logger)
        if substrate_report is None:
            return None
    else:
        logger.print(
            "[INFO] No substrate input detected. Substrate generation and substrate SDF copy will be skipped on both sides.")

    wt_report_dict = _run_mut_batch_side_workflow(
        cleaned_structure=wt_cleaned_structure,
        clean_report=wt_clean_report,
        input_msa=wt_input_msa,
        protein_name=wt_protein_name,
        msa_name=wt_msa_name,
        output_dir=wt_output_dir,
        substrate_names=resolved_substrate_names,
        substrate_report=substrate_report,
        substrate_dir=wt_output_dir if has_substrate else None,
        logger=logger,
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
    if wt_report_dict is None:
        return None

    mut_report_dict = _run_mut_batch_side_workflow(
        cleaned_structure=mut_cleaned_structure,
        clean_report=mut_clean_report,
        input_msa=mut_input_msa,
        protein_name=mut_protein_name,
        msa_name=mut_msa_name,
        output_dir=mut_output_dir,
        substrate_names=resolved_substrate_names,
        substrate_report=substrate_report,
        substrate_dir=mut_output_dir if has_substrate else None,
        logger=logger,
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
    if mut_report_dict is None:
        return None

    logger.print("[INFO] Mut_integrate workflow started")

    mut_integrate_report = integrate_mut_reports(
        mutclean_report=mutclean_report,
        wt_report_dict=wt_report_dict,
        mut_report_dict=mut_report_dict,
        strict=has_substrate,
        logger=logger,
    )
    if mut_integrate_report is None:
        return None

    logger.print("[INFO] Mut_batch workflow finished")


    return {
        "mut_clean_report": mutclean_report,
        "wt_report_dict": wt_report_dict,
        "mut_report_dict": mut_report_dict,
        "mut_integrate_report": mut_integrate_report,
    }