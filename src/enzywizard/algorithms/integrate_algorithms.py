from __future__ import annotations

from typing import Any, Dict, List, Tuple

from ..utils.logging_utils import Logger
from ..utils.integrate_utils import INTERACTION_ORDER,aa_class_statistics_to_list,aa_name_statistics_to_list,aa_ss_statistics_to_list,build_lookup_by_residue,get_clean_new_residue_list,get_disorder_membership_set,get_hydrophobic_cluster_membership_set,get_pocket_membership_set,interaction_one_hot,node_type_one_hot,residue_key,substrate_key
from ..utils.sequence_utils import normalize_aa_name_to_one_letter


def generate_integrate_report(overall_statistics: Dict[str, Any],integrated_graph: List[Dict[str, Any]]) -> Dict[str, Any]:
    return {
        "output_type": "enzywizard_integrate",
        "overall_statistics": overall_statistics,
        "integrated_graph": integrated_graph,
    }


def integrate_reports(report_dict: Dict[str, Dict[str, Any]],strict: bool,logger: Logger,) -> Dict[str, Any] | None:
    if not isinstance(report_dict, dict):
        logger.print("[ERROR] report_dict must be a dict.")
        return None

    clean_report = report_dict.get("enzywizard_clean")
    if not isinstance(clean_report, dict):
        logger.print("[ERROR] enzywizard_clean report is required.")
        return None

    if strict:
        required_types = [
            "enzywizard_clean",
            "enzywizard_aaprops",
            "enzywizard_hydrocluster",
            "enzywizard_energy",
            "enzywizard_flexibility",
            "enzywizard_disorder",
            "enzywizard_conservation",
            "enzywizard_embedding",
            "enzywizard_pocket",
            "enzywizard_substrate",
            "enzywizard_dock",
            "enzywizard_interaction",
        ]
        for output_type in required_types:
            if output_type not in report_dict:
                logger.print(f"[ERROR] Missing required report in strict mode: {output_type}")
                return None

    overall_statistics = build_overall_statistics(report_dict, strict, logger)
    if overall_statistics is None:
        return None

    integrated_graph = build_integrated_graph(report_dict, strict, logger)
    if integrated_graph is None:
        return None

    return generate_integrate_report(overall_statistics, integrated_graph)


def build_overall_statistics(report_dict: Dict[str, Dict[str, Any]],strict: bool,logger: Logger) -> Dict[str, Any] | None:
    overall_statistics: Dict[str, Any] = {}

    aaprops_report = report_dict.get("enzywizard_aaprops")
    hydro_report = report_dict.get("enzywizard_hydrocluster")
    disorder_report = report_dict.get("enzywizard_disorder")
    pocket_report = report_dict.get("enzywizard_pocket")
    energy_report = report_dict.get("enzywizard_energy")
    dock_report = report_dict.get("enzywizard_dock")
    interaction_report = report_dict.get("enzywizard_interaction")

    if (
        aaprops_report is None
        and hydro_report is None
        and disorder_report is None
        and pocket_report is None
        and energy_report is None
        and dock_report is None
        and interaction_report is None
    ):
        if strict:
            logger.print("[ERROR] overall_statistics cannot be empty in strict mode.")
            return None
        return overall_statistics

    if aaprops_report is not None:
        stats = aaprops_report.get("aa_props_statistics", {})
        aa_name_count = aa_name_statistics_to_list(stats.get("aa_name_statistics"), logger)
        if aa_name_count is None:
            return None
        aa_class_count = aa_class_statistics_to_list(stats.get("aa_class_statistics"), logger)
        if aa_class_count is None:
            return None
        aa_ss_count = aa_ss_statistics_to_list(stats.get("aa_ss_statistics"), logger)
        if aa_ss_count is None:
            return None

        overall_statistics["aa_name_count"] = aa_name_count
        overall_statistics["aa_class_count"] = aa_class_count
        overall_statistics["aa_ss_count"] = aa_ss_count
    elif strict:
        logger.print("[ERROR] enzywizard_aaprops is required in strict mode.")
        return None

    if hydro_report is not None:
        hydro_stats = hydro_report.get("hydrophobic_cluster_statistics", {})
        overall_statistics["hydrophobic_cluster_count"] = hydro_stats["cluster_num"]
        overall_statistics["max_hydrophobic_cluster_area"] = hydro_stats["max_cluster_area"]
        overall_statistics["total_hydrophobic_cluster_area"] = hydro_stats["total_cluster_area"]
    elif strict:
        logger.print("[ERROR] enzywizard_hydrocluster is required in strict mode.")
        return None

    if disorder_report is not None:
        disorder_stats = disorder_report.get("disorder_region_statistics", {})
        overall_statistics["disorder_region_count"] = disorder_stats["region_num"]
        overall_statistics["max_disorder_region_length"] = disorder_stats["max_region_length"]
        overall_statistics["total_disorder_region_length"] = disorder_stats["total_region_length"]
    elif strict:
        logger.print("[ERROR] enzywizard_disorder is required in strict mode.")
        return None

    if pocket_report is not None:
        pocket_stats = pocket_report.get("pocket_region_statistics", {})
        overall_statistics["pocket_region_count"] = pocket_stats["pocket_num"]
        overall_statistics["max_pocket_region_volume"] = pocket_stats["max_pocket_volume"]
        overall_statistics["total_pocket_region_volume"] = pocket_stats["total_pocket_volume"]
    elif strict:
        logger.print("[ERROR] enzywizard_pocket is required in strict mode.")
        return None

    if energy_report is not None:
        energy_terms = energy_report.get("energy_terms", {})
        overall_statistics["total_potential_energy"] = energy_terms["total_potential_energy"]
        overall_statistics["harmonic_bond_force"] = energy_terms["harmonic_bond_force"]
        overall_statistics["harmonic_angle_force"] = energy_terms["harmonic_angle_force"]
        overall_statistics["custom_bond_force"] = energy_terms["custom_bond_force"]
        overall_statistics["custom_torsion_force"] = energy_terms["custom_torsion_force"]
        overall_statistics["custom_nonbonded_force"] = energy_terms["custom_nonbonded_force"]
        overall_statistics["nonbonded_force"] = energy_terms["nonbonded_force"]
        overall_statistics["periodic_torsion_force"] = energy_terms["periodic_torsion_force"]
        overall_statistics["cmap_torsion_force"] = energy_terms["cmap_torsion_force"]
    elif strict:
        logger.print("[ERROR] enzywizard_energy is required in strict mode.")
        return None

    if dock_report is not None:
        overall_statistics["docking_score"] = dock_report["docked_result"]["docking_score"]
    elif strict:
        logger.print("[ERROR] enzywizard_dock is required in strict mode.")
        return None

    if interaction_report is not None:
        count_block = interaction_report["interactions_statistics"]["overall"]["count"]
        overall_statistics["hbond_count"] = count_block["HBOND"]
        overall_statistics["ionic_count"] = count_block["IONIC"]
        overall_statistics["vdw_count"] = count_block["VDW"]
        overall_statistics["pipistack_count"] = count_block["PIPISTACK"]
        overall_statistics["pication_count"] = count_block["PICATION"]
        overall_statistics["ssbond_count"] = count_block["SSBOND"]
    elif strict:
        logger.print("[ERROR] enzywizard_interaction is required in strict mode.")
        return None

    return overall_statistics


def build_integrated_graph(report_dict: Dict[str, Dict[str, Any]],strict: bool,logger: Logger) -> List[Dict[str, Any]] | None:
    interaction_report = report_dict.get("enzywizard_interaction")
    if interaction_report is None:
        if strict:
            logger.print("[ERROR] enzywizard_interaction is required in strict mode.")
            return None
        return []

    interactions = interaction_report.get("interactions")
    if not isinstance(interactions, list):
        logger.print("[ERROR] interactions must be a list.")
        return None

    if len(interactions) == 0:
        return []

    amino_acid_node_list = build_amino_acid_nodes(report_dict, strict, logger)
    if amino_acid_node_list is None:
        return None

    substrate_node_list = build_substrate_nodes(report_dict, strict, logger)
    if substrate_node_list is None:
        return None

    all_node_list = amino_acid_node_list + substrate_node_list

    for i, node in enumerate(all_node_list, start=1):
        node["node_id"] = i

    node_entry_list: List[Dict[str, Any]] = []
    for node in all_node_list:
        node_entry_list.append({
            "node_1": node
        })

    node_lookup = build_integrated_node_lookup(all_node_list, logger)
    if node_lookup is None:
        return None

    edge_entry_list = build_edge_entries_from_interactions(interactions, node_lookup, strict, logger)
    if edge_entry_list is None:
        return None

    return node_entry_list + edge_entry_list


def build_amino_acid_nodes(report_dict: Dict[str, Dict[str, Any]],strict: bool,logger: Logger,) -> List[Dict[str, Any]] | None:
    clean_report = report_dict.get("enzywizard_clean")
    if clean_report is None:
        logger.print("[ERROR] Missing clean report.")
        return None

    aaprops_report = report_dict.get("enzywizard_aaprops")
    hydro_report = report_dict.get("enzywizard_hydrocluster")
    flexibility_report = report_dict.get("enzywizard_flexibility")
    disorder_report = report_dict.get("enzywizard_disorder")
    conservation_report = report_dict.get("enzywizard_conservation")
    embedding_report = report_dict.get("enzywizard_embedding")
    pocket_report = report_dict.get("enzywizard_pocket")

    clean_residue_list = get_clean_new_residue_list(clean_report, logger)
    if clean_residue_list is None:
        return None

    aaprops_lookup = None
    if aaprops_report is not None:
        aaprops_lookup = build_lookup_by_residue(aaprops_report["aa_props"], "aa_id", "aa_name", logger)
        if aaprops_lookup is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing aaprops report in strict mode.")
        return None

    flexibility_lookup = None
    if flexibility_report is not None:
        flexibility_lookup = build_lookup_by_residue(flexibility_report["protein_rmsf"], "aa_id", "aa_name", logger)
        if flexibility_lookup is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing flexibility report in strict mode.")
        return None

    conservation_lookup = None
    if conservation_report is not None:
        conservation_lookup = build_lookup_by_residue(conservation_report["conservation_scores"], "aa_id", "aa_name", logger)
        if conservation_lookup is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing conservation report in strict mode.")
        return None

    embedding_lookup = None
    if embedding_report is not None:
        embedding_lookup = build_lookup_by_residue(embedding_report["embeddings"], "aa_id", "aa_name", logger)
        if embedding_lookup is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing embedding report in strict mode.")
        return None

    hydrophobic_membership = None
    if hydro_report is not None:
        hydrophobic_membership = get_hydrophobic_cluster_membership_set(hydro_report, logger)
        if hydrophobic_membership is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing hydrocluster report in strict mode.")
        return None

    disorder_membership = None
    if disorder_report is not None:
        disorder_membership = get_disorder_membership_set(disorder_report, logger)
        if disorder_membership is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing disorder report in strict mode.")
        return None

    pocket_membership = None
    if pocket_report is not None:
        pocket_membership = get_pocket_membership_set(pocket_report, logger)
        if pocket_membership is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing pocket report in strict mode.")
        return None

    node_list: List[Dict[str, Any]] = []

    for residue in clean_residue_list:
        aa_index = residue["aa_id"]
        aa_name = normalize_aa_name_to_one_letter(residue["aa_name"])
        key = residue_key(aa_index, aa_name)

        drop_node = False

        node: Dict[str, Any] = {}
        node["node_id"] = 0
        node["node_type"] = "amino_acid"
        node["node_type_one_hot"] = node_type_one_hot("amino_acid")
        node["aa_index"] = aa_index
        node["aa_name"] = aa_name

        if aaprops_lookup is not None:
            aa_item = aaprops_lookup.get(key)
            if aa_item is None:
                if strict:
                    logger.print(f"[ERROR] Missing aa_props entry for residue {aa_index} {aa_name}.")
                    return None
                logger.print(f"[WARNING] Missing aa_props entry for residue {aa_index} {aa_name}. Node dropped.")
                drop_node = True
            else:
                node["aa_coord"] = aa_item["aa_coord"]
                node["aa_class"] = aa_item["aa_class"]
                node["aa_ss"] = aa_item["aa_ss"]
                node["aa_rsa"] = aa_item["aa_rsa"]
                node["aa_phi"] = aa_item["aa_phi"]
                node["aa_psi"] = aa_item["aa_psi"]
                node["aa_net_charge"] = aa_item["aa_net_charge"]
                node["aa_pka"] = aa_item["aa_pka"]
                node["aa_volume"] = aa_item["aa_volume"]
                node["aa_hydrophobicity"] = aa_item["aa_hydrophobicity"]
                node["aa_molecular_weight"] = aa_item["aa_molecular_weight"]
                node["aa_pi"] = aa_item["aa_pi"]
                node["aa_name_one_hot"] = aa_item["aa_name_one_hot"]
                node["aa_class_one_hot"] = aa_item["aa_class_one_hot"]
                node["aa_ss_one_hot"] = aa_item["aa_ss_one_hot"]

        if drop_node:
            continue

        if hydrophobic_membership is not None:
            node["is_in_hydrophobic_cluster"] = (key in hydrophobic_membership)

        if flexibility_lookup is not None:
            rmsf_item = flexibility_lookup.get(key)
            if rmsf_item is None:
                if strict:
                    logger.print(f"[ERROR] Missing RMSF entry for residue {aa_index} {aa_name}.")
                    return None
                logger.print(f"[WARNING] Missing RMSF entry for residue {aa_index} {aa_name}. Node dropped.")
                continue
            node["rmsf"] = rmsf_item["rmsf"]

        if disorder_membership is not None:
            node["is_in_disorder_region"] = (key in disorder_membership)

        if conservation_lookup is not None:
            cons_item = conservation_lookup.get(key)
            if cons_item is None:
                if strict:
                    logger.print(f"[ERROR] Missing conservation entry for residue {aa_index} {aa_name}.")
                    return None
                logger.print(f"[WARNING] Missing conservation entry for residue {aa_index} {aa_name}. Node dropped.")
                continue
            node["conservation_score"] = cons_item["conservation_score"]

        if embedding_lookup is not None:
            emb_item = embedding_lookup.get(key)
            if emb_item is None:
                if strict:
                    logger.print(f"[ERROR] Missing embedding entry for residue {aa_index} {aa_name}.")
                    return None
                logger.print(f"[WARNING] Missing embedding entry for residue {aa_index} {aa_name}. Node dropped.")
                continue
            node["embedding"] = emb_item["embedding"]

        if pocket_membership is not None:
            node["is_in_pocket"] = (key in pocket_membership)

        node = reorder_amino_acid_node(node)
        node_list.append(node)

    return node_list


def build_substrate_nodes(report_dict: Dict[str, Dict[str, Any]],strict: bool,logger: Logger) -> List[Dict[str, Any]] | None:
    substrate_report = report_dict.get("enzywizard_substrate")
    dock_report = report_dict.get("enzywizard_dock")

    if substrate_report is None or dock_report is None:
        if strict:
            logger.print("[ERROR] substrate and dock reports are required in strict mode.")
            return None
        return []

    substrate_lookup: Dict[str, Dict[str, Any]] = {}
    for item in substrate_report["substrates"]:
        substrate_name = item["substrate_name"]
        if substrate_name in substrate_lookup:
            logger.print(f"[ERROR] Duplicate substrate_name in substrate report: {substrate_name}")
            return None
        substrate_lookup[substrate_name] = item

    docked_substrate_list = dock_report["docked_result"]["docked_substrates"]
    node_list: List[Dict[str, Any]] = []

    for substrate_index, dock_item in enumerate(docked_substrate_list, start=1):
        substrate_name = dock_item["substrate_name"]

        substrate_item = substrate_lookup.get(substrate_name)
        if substrate_item is None:
            if strict:
                logger.print(f"[ERROR] Missing substrate entry for docked substrate: {substrate_name}")
                return None
            logger.print(f"[WARNING] Missing substrate entry for docked substrate: {substrate_name}. Node dropped.")
            continue

        node: Dict[str, Any] = {}
        node["node_id"] = 0
        node["node_type"] = "substrate"
        node["node_type_one_hot"] = node_type_one_hot("substrate")
        node["substrate_index"] = substrate_index
        node["substrate_name"] = substrate_name
        node["smiles"] = substrate_item["smiles"]
        node["num_atoms"] = substrate_item["num_atoms"]
        node["mol_weight"] = substrate_item["mol_weight"]
        node["logp"] = substrate_item["logp"]
        node["docked_center_coord"] = dock_item["docked_center_coord"]
        node["fingerprint"] = substrate_item["fingerprint"]

        node = reorder_substrate_node(node)
        node_list.append(node)

    return node_list


def build_integrated_node_lookup(all_node_list: List[Dict[str, Any]],logger: Logger,) -> Dict[Tuple[Any, ...], Dict[str, Any]] | None:
    lookup: Dict[Tuple[Any, ...], Dict[str, Any]] = {}

    for node in all_node_list:
        node_type = node.get("node_type")
        if node_type == "amino_acid":
            key = ("amino_acid", node["aa_index"], node["aa_name"])
        elif node_type == "substrate":
            key = ("substrate", node["substrate_name"])
        else:
            logger.print("[ERROR] Unknown node_type while building node lookup.")
            return None

        if key in lookup:
            logger.print(f"[ERROR] Duplicate integrated node key: {key}")
            return None

        lookup[key] = node

    return lookup


def build_edge_entries_from_interactions(interactions: List[Dict[str, Any]],node_lookup: Dict[Tuple[str, int, str], Dict[str, Any]],strict: bool,logger: Logger) -> List[Dict[str, Any]] | None:
    merged_edge_count: Dict[Tuple[str, int, int], int] = {}
    merged_edge_nodes: Dict[Tuple[str, int, int], Tuple[Dict[str, Any], Dict[str, Any]]] = {}

    for item in interactions:
        interaction = item["interaction"]
        node1 = item["node1"]
        node2 = item["node2"]

        integrated_node_1 = resolve_interaction_node_to_integrated_node(node1, node_lookup, logger)
        if integrated_node_1 is None:
            if strict:
                logger.print("[ERROR] Failed to resolve interaction node1 in strict mode.")
                return None
            logger.print("[WARNING] Interaction node1 cannot be resolved. Edge skipped.")
            continue

        integrated_node_2 = resolve_interaction_node_to_integrated_node(node2, node_lookup, logger)
        if integrated_node_2 is None:
            if strict:
                logger.print("[ERROR] Failed to resolve interaction node2 in strict mode.")
                return None
            logger.print("[WARNING] Interaction node2 cannot be resolved. Edge skipped.")
            continue

        node_id_1 = integrated_node_1["node_id"]
        node_id_2 = integrated_node_2["node_id"]

        if node_id_1 <= node_id_2:
            left_node = integrated_node_1
            right_node = integrated_node_2
            left_id = node_id_1
            right_id = node_id_2
        else:
            left_node = integrated_node_2
            right_node = integrated_node_1
            left_id = node_id_2
            right_id = node_id_1

        merge_key = (interaction, left_id, right_id)

        merged_edge_count[merge_key] = merged_edge_count.get(merge_key, 0) + 1
        merged_edge_nodes[merge_key] = (left_node, right_node)

    edge_entry_list: List[Dict[str, Any]] = []

    for interaction, node_id_1, node_id_2 in sorted(merged_edge_count.keys(), key=lambda x: (x[1], x[2], INTERACTION_ORDER.index(x[0]))):
        node_1, node_2 = merged_edge_nodes[(interaction, node_id_1, node_id_2)]

        edge = {
            "interaction": interaction,
            "interaction_one_hot": interaction_one_hot(interaction),
            "interaction_count": merged_edge_count[(interaction, node_id_1, node_id_2)],
        }

        edge_entry_list.append({
            "edge": edge,
            "node_1": node_1,
            "node_2": node_2,
        })

    return edge_entry_list


def resolve_interaction_node_to_integrated_node(interaction_node: Dict[str, Any],node_lookup: Dict[Tuple[Any, ...], Dict[str, Any]],logger: Logger,) -> Dict[str, Any] | None:
    node_type = interaction_node.get("node_type")

    if node_type == "amino_acid":
        aa_index = interaction_node.get("aa_index")
        aa_name = interaction_node.get("aa_name")

        if not isinstance(aa_index, int) or not isinstance(aa_name, str):
            logger.print("[ERROR] Invalid amino acid interaction node.")
            return None

        aa_name_norm = normalize_aa_name_to_one_letter(aa_name)

        return node_lookup.get(("amino_acid", aa_index, aa_name_norm))

    if node_type == "substrate":
        substrate_name = interaction_node.get("substrate_name")
        if not isinstance(substrate_name, str) or substrate_name.strip() == "":
            logger.print("[ERROR] Invalid substrate interaction node.")
            return None

        substrate_name_clean = substrate_name.strip()

        node = node_lookup.get(("substrate", substrate_name_clean))
        if node is not None:
            return node

        if substrate_name_clean.startswith("docked_"):
            substrate_name_clean = substrate_name_clean[len("docked_"):]

        node = node_lookup.get(("substrate", substrate_name_clean))
        if node is not None:
            return node

        logger.print(f"[WARNING] Failed to match substrate node by substrate_name: {substrate_name}")
        return None


def reorder_amino_acid_node(node: Dict[str, Any]) -> Dict[str, Any]:
    ordered: Dict[str, Any] = {}

    ordered["node_id"] = node["node_id"]
    ordered["node_type"] = node["node_type"]
    ordered["node_type_one_hot"] = node["node_type_one_hot"]

    ordered["aa_index"] = node["aa_index"]
    ordered["aa_name"] = node["aa_name"]

    if "aa_name_one_hot" in node:
        ordered["aa_name_one_hot"] = node["aa_name_one_hot"]

    if "aa_coord" in node:
        ordered["aa_coord"] = node["aa_coord"]

    if "aa_class" in node:
        ordered["aa_class"] = node["aa_class"]

    if "aa_class_one_hot" in node:
        ordered["aa_class_one_hot"] = node["aa_class_one_hot"]

    if "aa_ss" in node:
        ordered["aa_ss"] = node["aa_ss"]

    if "aa_ss_one_hot" in node:
        ordered["aa_ss_one_hot"] = node["aa_ss_one_hot"]

    if "aa_rsa" in node:
        ordered["aa_rsa"] = node["aa_rsa"]
    if "aa_phi" in node:
        ordered["aa_phi"] = node["aa_phi"]
    if "aa_psi" in node:
        ordered["aa_psi"] = node["aa_psi"]
    if "aa_net_charge" in node:
        ordered["aa_net_charge"] = node["aa_net_charge"]
    if "aa_pka" in node:
        ordered["aa_pka"] = node["aa_pka"]
    if "aa_volume" in node:
        ordered["aa_volume"] = node["aa_volume"]
    if "aa_hydrophobicity" in node:
        ordered["aa_hydrophobicity"] = node["aa_hydrophobicity"]
    if "aa_molecular_weight" in node:
        ordered["aa_molecular_weight"] = node["aa_molecular_weight"]
    if "aa_pi" in node:
        ordered["aa_pi"] = node["aa_pi"]

    if "rmsf" in node:
        ordered["rmsf"] = node["rmsf"]
    if "conservation_score" in node:
        ordered["conservation_score"] = node["conservation_score"]

    if "embedding" in node:
        ordered["embedding"] = node["embedding"]

    if "is_in_hydrophobic_cluster" in node:
        ordered["is_in_hydrophobic_cluster"] = node["is_in_hydrophobic_cluster"]
    if "is_in_disorder_region" in node:
        ordered["is_in_disorder_region"] = node["is_in_disorder_region"]
    if "is_in_pocket" in node:
        ordered["is_in_pocket"] = node["is_in_pocket"]

    return ordered


def reorder_substrate_node(node: Dict[str, Any]) -> Dict[str, Any]:
    ordered: Dict[str, Any] = {}

    ordered["node_id"] = node["node_id"]
    ordered["node_type"] = node["node_type"]
    ordered["node_type_one_hot"] = node["node_type_one_hot"]

    ordered["substrate_index"] = node["substrate_index"]
    ordered["substrate_name"] = node["substrate_name"]

    if "smiles" in node:
        ordered["smiles"] = node["smiles"]
    if "num_atoms" in node:
        ordered["num_atoms"] = node["num_atoms"]
    if "mol_weight" in node:
        ordered["mol_weight"] = node["mol_weight"]
    if "logp" in node:
        ordered["logp"] = node["logp"]
    if "docked_center_coord" in node:
        ordered["docked_center_coord"] = node["docked_center_coord"]
    if "fingerprint" in node:
        ordered["fingerprint"] = node["fingerprint"]

    return ordered