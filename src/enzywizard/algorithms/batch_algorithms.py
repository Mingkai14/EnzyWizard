from __future__ import annotations

from pathlib import Path
from typing import Dict, Any
from openmm.app import Modeller

from ..utils.logging_utils import Logger
from ..utils.batch_utils import build_batch_output_paths

from ..utils.IO_utils import structure_to_pdbfile,load_dssp,load_msa,write_msa,write_hmm,save_substrate_structures,load_sdf_mol_3d

from ..utils.sequence_utils import check_msa, clean_msa_to_sto
from ..utils.structure_utils import structure_has_hydrogen, get_fasta_dict_from_structure
from ..utils.interaction_utils import filter_valid_docked_substrates

from ..algorithms.clean_algorithms import clean_structure_to_single_chain_A,generate_clean_report,check_cleaned_structure,validate_clean_mapping_coordinates

from ..algorithms.aaprops_algorithms import calculate_aa_props,calculate_aa_props_statistics,generate_aaprops_report

from ..algorithms.hydrocluster_algorithms import compute_hydrophobic_clusters,generate_hydrocluster_report

from ..algorithms.energy_algorithms import compute_energy_terms,generate_energy_report

from ..algorithms.flexibility_algorithms import compute_protein_rmsf,generate_flexibility_report

from ..algorithms.disorder_algorithms import compute_disordered_regions,generate_disorder_report

from ..algorithms.conservation_algorithms import compute_conservation_scores,generate_conservation_report

from ..algorithms.embedding_algorithms import generate_embedding,generate_embedding_report

from ..algorithms.pocket_algorithms import compute_pockets,generate_pocket_report

from ..algorithms.substrate_algorithms import get_substrate_dict_list_from_input,get_completed_smiles_list,get_substrate_feature_list,generate_substrate_report

from ..algorithms.dock_algorithms import dock_multiple_substrates_from_structure,save_docking_results_and_generate_dock_report

from ..algorithms.interaction_algorithms import calculate_all_interaction_network,summarize_interaction_counts,generate_interaction_report

from ..algorithms.integrate_algorithms import integrate_reports

from ..utils.substrate_utils import build_docked_mol_from_atom_info

from ..utils.IO_utils import load_protein_structure

def run_batch_workflow(
    cleaned_input_path: str | Path,
    input_msa: str | Path,
    substrate_names: str | None,
    protein_name: str,
    msa_name: str,
    output_dir: str | Path,
    logger: Logger,
    cutoff_area: float = 10.0,
    minimize_energy: bool = True,
    minimization_iteration: int = 2000,
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
    cleaned_input_path = Path(cleaned_input_path)
    input_msa = Path(input_msa)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    path_dict = build_batch_output_paths(protein_name=protein_name, msa_name=msa_name, output_dir=output_dir)

    report_dict: Dict[str, Dict[str, Any]] = {}

    has_substrate = isinstance(substrate_names, str) and substrate_names.strip() != ""

    logger.print("[INFO] Batch workflow started from cleaned input structure")

    try:
        original_structure = load_protein_structure(cleaned_input_path, protein_name, logger)
    except Exception:
        original_structure = None

    if original_structure is None:
        logger.print(f"[ERROR] Failed to load structure: {cleaned_input_path}")
        return None

    if not check_cleaned_structure(original_structure, logger):
        logger.print("[ERROR] Input structure is not a valid cleaned structure.")
        return None

    if not structure_has_hydrogen(original_structure, logger):
        logger.print("[ERROR] Input cleaned structure does not contain hydrogen atoms.")
        return None

    clean_result = clean_structure_to_single_chain_A(original_structure, logger)
    if clean_result is None:
        return None

    cleaned_structure, mapping_old_to_new, clean_stats = clean_result


    if not validate_clean_mapping_coordinates(original_structure, cleaned_structure, mapping_old_to_new, logger):
        return None

    clean_report = generate_clean_report(
        original_structure,
        cleaned_structure,
        mapping_old_to_new,
        clean_stats,
        logger,
    )
    if clean_report is None:
        return None
    report_dict["enzywizard_clean"] = clean_report

    cleaned_pdbfile = structure_to_pdbfile(cleaned_structure, logger, protein_name=protein_name)
    if cleaned_pdbfile is None:
        return None

    try:
        cleaned_modeller = Modeller(cleaned_pdbfile.topology, cleaned_pdbfile.positions)
    except Exception as e:
        logger.print(f"[ERROR] Failed to build OpenMM Modeller: {e}")
        return None
    logger.print("[INFO] Cleaned OpenMM Modeller loaded")

    sequence_dict = get_fasta_dict_from_structure(cleaned_structure, logger, header=protein_name)
    if sequence_dict is None:
        return None
    logger.print("[INFO] Cleaned sequence prepared")

    logger.print("[INFO] Aaprops calculation started")
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

    logger.print("[INFO] Hydrocluster calculation started")
    clusters = compute_hydrophobic_clusters(cleaned_structure, logger, cutoff_area=cutoff_area)
    if clusters is None:
        return None

    hydrocluster_report = generate_hydrocluster_report(clusters, cleaned_structure, logger)
    if hydrocluster_report is None:
        return None
    report_dict["enzywizard_hydrocluster"] = hydrocluster_report

    logger.print("[INFO] Energy calculation started")
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

    logger.print("[INFO] Flexibility calculation started")
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

    logger.print("[INFO] Disorder calculation started")
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

    logger.print("[INFO] Conservation calculation started")
    msa_list = load_msa(input_msa, logger)
    if msa_list is None:
        return None

    if not check_msa(input_msa, sequence_dict, msa_list, logger):
        return None

    cleaned_msa_list = clean_msa_to_sto(msa_list, logger)
    if cleaned_msa_list is None:
        return None

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

    logger.print("[INFO] Embedding calculation started")
    embeddings = generate_embedding(sequence_dict, logger, model_name=embedding_model_name)
    if embeddings is None:
        return None

    embedding_report = generate_embedding_report(embeddings)
    report_dict["enzywizard_embedding"] = embedding_report

    logger.print("[INFO] Pocket calculation started")
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

    ligand_mol_list = []
    substrate_name_list = []

    if has_substrate:
        logger.print("[INFO] Substrate calculation started")
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

        if not save_substrate_structures(substrate_feature_list, output_dir, logger):
            return None
        logger.print(f"[INFO] Substrate structures saved: {output_dir}")

        substrate_report = generate_substrate_report(substrate_feature_list, logger)
        if substrate_report is None:
            return None
        report_dict["enzywizard_substrate"] = substrate_report


        logger.print("[INFO] Docking workflow started")
        docking_result_list = dock_multiple_substrates_from_structure(
            struct=cleaned_structure,
            substrate_names=resolved_substrate_names,
            substrate_dir=output_dir,
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

        logger.print(f"[INFO] Loaded {len(ligand_mol_list)} docked substrate Mol(3D) object(s)")
    else:
        logger.print("[INFO] No substrate input detected. Skipping substrate and docking workflows.")

    logger.print("[INFO] Interaction workflow started")

    if has_substrate:
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
        logger.print(f"[INFO] Valid docked substrate count: {len(valid_ligand_mol_list)}")
    else:
        valid_substrate_name_list = []
        valid_ligand_mol_list = []
        logger.print("[INFO] No substrate input detected. Only intra-protein interactions will be calculated.")

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
        logger.print("[ERROR] Failed to calculate interaction network.")
        return None

    interaction_statistics = summarize_interaction_counts(interaction_list=interaction_list, logger=logger)
    if interaction_statistics is None:
        return None

    interaction_report = generate_interaction_report(
        interaction_list=interaction_list,
        interaction_statistics=interaction_statistics,
    )
    report_dict["enzywizard_interaction"] = interaction_report

    logger.print("[INFO] Integrate workflow started")
    integrate_strict = has_substrate
    integrate_report = integrate_reports(report_dict, integrate_strict, logger)
    if integrate_report is None:
        return None

    integrated_graph = integrate_report.get("integrated_graph")
    if not isinstance(integrated_graph, list):
        logger.print("[ERROR] integrated_graph missing in integrate report.")
        return None

    return {
        "integrate_report": integrate_report,
        "report_dict": report_dict,
    }