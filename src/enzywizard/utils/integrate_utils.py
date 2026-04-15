from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple
import json
import re

from ..utils.logging_utils import Logger
from ..utils.IO_utils import write_json_from_dict_inline_leaf_lists
from ..utils.common_utils import get_optimized_filename
from ..resources.aa_resources import AA_20NAME_INDEX, AA_8CLASSES, DSSP_8STATE_INDEX
from ..utils.sequence_utils import normalize_aa_name_to_one_letter


SUPPORTED_OUTPUT_TYPES = {
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
}

INTERACTION_ORDER = ["HBOND", "IONIC", "VDW", "PIPISTACK", "PICATION", "SSBOND"]
NODE_TYPE_ORDER = ["amino_acid", "substrate"]
AA_CLASS_ORDER = [item[0] for item in AA_8CLASSES]
AA_NAME_ORDER = [aa for aa, _ in sorted(AA_20NAME_INDEX.items(), key=lambda x: x[1])]
AA_SS_ORDER = [ss for ss, _ in sorted(DSSP_8STATE_INDEX.items(), key=lambda x: x[1])]


def load_json_file(path: str | Path, logger: Logger) -> Dict[str, Any] | None:
    try:
        p = Path(path)
        if not p.exists() or not p.is_file():
            logger.print(f"[ERROR] JSON file not found: {p}")
            return None

        with p.open("r", encoding="utf-8") as f:
            data = json.load(f)

        if not isinstance(data, dict):
            logger.print(f"[ERROR] JSON top-level object must be a dict: {p}")
            return None

        return data

    except Exception as e:
        logger.print(f"[ERROR] Failed to load JSON file {path}: {e}")
        return None


def list_json_files(input_dir: str | Path, logger: Logger) -> List[Path] | None:
    try:
        p = Path(input_dir)
        if not p.exists() or not p.is_dir():
            logger.print(f"[ERROR] Invalid input_dir: {p}")
            return None

        json_files = sorted([x for x in p.iterdir() if x.is_file() and x.suffix.lower() == ".json"])
        return json_files

    except Exception as e:
        logger.print(f"[ERROR] Failed to list JSON files from {input_dir}: {e}")
        return None

def find_unique_clean_report_path(json_path_list: List[Path], logger: Logger) -> Path | None:
    if not isinstance(json_path_list, list):
        logger.print("[ERROR] json_path_list must be a list.")
        return None

    clean_report_path_list: List[Path] = []

    for json_path in json_path_list:
        data = load_json_file(json_path, logger)
        if data is None:
            return None

        output_type = get_supported_output_type(data, logger)
        if output_type is None:
            return None

        if output_type == "enzywizard_clean":
            clean_report_path_list.append(json_path)

    if len(clean_report_path_list) == 0:
        logger.print("[ERROR] No enzywizard_clean report found in input_dir.")
        return None

    if len(clean_report_path_list) > 1:
        logger.print(f"[ERROR] Multiple enzywizard_clean reports found in input_dir: {len(clean_report_path_list)}")
        return None

    return clean_report_path_list[0]


def extract_protein_name_from_clean_report_path(clean_report_path: str | Path, logger: Logger) -> str | None:
    try:
        name = Path(clean_report_path).name
        m = re.fullmatch(r"clean_report_(.+)\.json", name)
        if m is None:
            logger.print(f"[ERROR] clean_report file name must match clean_report_{{protein_name}}.json: {name}")
            return None

        protein_name = m.group(1).strip()
        if protein_name == "":
            logger.print(f"[ERROR] Invalid protein name in clean_report file name: {name}")
            return None

        return protein_name

    except Exception as e:
        logger.print(f"[ERROR] Failed to parse clean_report file name: {e}")
        return None


def save_integrate_json(report: Dict[str, Any] | List[Any], output_path: str | Path, logger: Logger) -> bool:
    try:
        write_json_from_dict_inline_leaf_lists(report, output_path)
        return True
    except Exception as e:
        logger.print(f"[ERROR] Failed to save integrate JSON: {e}")
        return False

def split_integrated_graph_entries(
    integrated_graph: List[Dict[str, Any]],
    logger: Logger,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]] | None:
    if not isinstance(integrated_graph, list):
        logger.print("[ERROR] integrated_graph must be a list.")
        return None

    node_list: List[Dict[str, Any]] = []
    edge_list: List[Dict[str, Any]] = []

    for item in integrated_graph:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid integrated_graph item.")
            return None

        if "node_1" in item and "edge" not in item and "node_2" not in item:
            node_list.append(item["node_1"])
        elif "edge" in item and "node_1" in item and "node_2" in item:
            edge_list.append(item)
        else:
            logger.print("[ERROR] Invalid integrated_graph entry structure.")
            return None

    return node_list, edge_list



def interaction_one_hot(interaction: str) -> List[int]:
    vec = [0] * len(INTERACTION_ORDER)
    if interaction in INTERACTION_ORDER:
        vec[INTERACTION_ORDER.index(interaction)] = 1
    return vec


def node_type_one_hot(node_type: str) -> List[int]:
    vec = [0] * len(NODE_TYPE_ORDER)
    if node_type in NODE_TYPE_ORDER:
        vec[NODE_TYPE_ORDER.index(node_type)] = 1
    return vec


def aa_name_statistics_to_list(statistics: Dict[str, int], logger: Logger) -> List[int] | None:
    if not isinstance(statistics, dict):
        logger.print("[ERROR] aa_name_statistics must be a dict.")
        return None

    out = [0] * len(AA_NAME_ORDER)
    for aa in AA_NAME_ORDER:
        value = statistics.get(aa)
        if not isinstance(value, int):
            logger.print(f"[ERROR] Missing or invalid aa_name_statistics for {aa}.")
            return None
        out[AA_20NAME_INDEX[aa]] = value
    return out


def aa_class_statistics_to_list(statistics: Dict[str, int], logger: Logger) -> List[int] | None:
    if not isinstance(statistics, dict):
        logger.print("[ERROR] aa_class_statistics must be a dict.")
        return None

    out = [0] * len(AA_CLASS_ORDER)
    for i, cls in enumerate(AA_CLASS_ORDER):
        value = statistics.get(cls)
        if not isinstance(value, int):
            logger.print(f"[ERROR] Missing or invalid aa_class_statistics for {cls}.")
            return None
        out[i] = value
    return out


def aa_ss_statistics_to_list(statistics: Dict[str, int], logger: Logger) -> List[int] | None:
    if not isinstance(statistics, dict):
        logger.print("[ERROR] aa_ss_statistics must be a dict.")
        return None

    out = [0] * len(AA_SS_ORDER)
    for ss in AA_SS_ORDER:
        value = statistics.get(ss)
        if not isinstance(value, int):
            logger.print(f"[ERROR] Missing or invalid aa_ss_statistics for {ss}.")
            return None
        out[DSSP_8STATE_INDEX[ss]] = value
    return out



def residue_key(aa_id: int, aa_name: str) -> Tuple[int, str]:
    aa_name_norm = normalize_aa_name_to_one_letter(aa_name)
    return int(aa_id), aa_name_norm


def residue_key_from_dict(data: Dict[str, Any], id_key: str, name_key: str) -> Tuple[int, str] | None:
    aa_id = data.get(id_key)
    aa_name = data.get(name_key)

    if not isinstance(aa_id, int):
        return None
    if not isinstance(aa_name, str) or aa_name.strip() == "":
        return None

    return residue_key(aa_id, aa_name)


def substrate_key(substrate_index: int, substrate_name: str) -> Tuple[int, str]:
    return int(substrate_index), str(substrate_name).strip()


def build_residue_membership_set(residue_list: List[Dict[str, Any]], logger: Logger) -> set[Tuple[int, str]] | None:
    if not isinstance(residue_list, list):
        logger.print("[ERROR] residue_list must be a list.")
        return None

    out: set[Tuple[int, str]] = set()

    for item in residue_list:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid residue item.")
            return None

        key = residue_key_from_dict(item, "aa_id", "aa_name")
        if key is None:
            logger.print("[ERROR] Invalid residue key in residue list.")
            return None
        out.add(key)

    return out


def build_lookup_by_residue(
    item_list: List[Dict[str, Any]],
    id_key: str,
    name_key: str,
    logger: Logger
) -> Dict[Tuple[int, str], Dict[str, Any]] | None:
    if not isinstance(item_list, list):
        logger.print("[ERROR] item_list must be a list.")
        return None

    lookup: Dict[Tuple[int, str], Dict[str, Any]] = {}

    for item in item_list:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid item in lookup list.")
            return None

        key = residue_key_from_dict(item, id_key, name_key)
        if key is None:
            logger.print("[ERROR] Invalid residue identifier in lookup list.")
            return None

        if key in lookup:
            logger.print(f"[ERROR] Duplicate residue entry found: {key}")
            return None

        lookup[key] = item

    return lookup


def get_hydrophobic_cluster_membership_set(hydro_report: Dict[str, Any], logger: Logger) -> set[Tuple[int, str]] | None:
    clusters = hydro_report.get("hydrophobic_cluster")
    if not isinstance(clusters, list):
        logger.print("[ERROR] Invalid hydrophobic_cluster.")
        return None

    out: set[Tuple[int, str]] = set()

    for cluster in clusters:
        if not isinstance(cluster, dict):
            logger.print("[ERROR] Invalid hydrophobic cluster entry.")
            return None
        residues = cluster.get("residues")
        membership = build_residue_membership_set(residues, logger)
        if membership is None:
            return None
        out |= membership

    return out


def get_disorder_membership_set(disorder_report: Dict[str, Any], logger: Logger) -> set[Tuple[int, str]] | None:
    regions = disorder_report.get("disorder_regions")
    if not isinstance(regions, list):
        logger.print("[ERROR] Invalid disorder_regions.")
        return None

    out: set[Tuple[int, str]] = set()

    for region in regions:
        if not isinstance(region, dict):
            logger.print("[ERROR] Invalid disorder region entry.")
            return None
        residues = region.get("residues")
        membership = build_residue_membership_set(residues, logger)
        if membership is None:
            return None
        out |= membership

    return out


def get_pocket_membership_set(pocket_report: Dict[str, Any], logger: Logger) -> set[Tuple[int, str]] | None:
    regions = pocket_report.get("pocket_regions")
    if not isinstance(regions, list):
        logger.print("[ERROR] Invalid pocket_regions.")
        return None

    out: set[Tuple[int, str]] = set()

    for region in regions:
        if not isinstance(region, dict):
            logger.print("[ERROR] Invalid pocket region entry.")
            return None
        residues = region.get("residues")
        membership = build_residue_membership_set(residues, logger)
        if membership is None:
            return None
        out |= membership

    return out


def get_clean_new_residue_list(clean_report: Dict[str, Any], logger: Logger) -> List[Dict[str, Any]] | None:
    mapping = clean_report.get("amino_acid_mapping_old_to_new")
    if not isinstance(mapping, list):
        logger.print("[ERROR] Invalid amino_acid_mapping_old_to_new in clean report.")
        return None

    if len(mapping) == 0:
        logger.print("[ERROR] amino_acid_mapping_old_to_new is empty.")
        return None

    new_residue_list: List[Dict[str, Any]] = []

    for item in mapping:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid mapping entry in clean report.")
            return None

        new_residue = item.get("new_residue")
        if not isinstance(new_residue, dict):
            logger.print("[ERROR] Missing new_residue in clean report.")
            return None

        aa_id = new_residue.get("aa_id")
        aa_name = new_residue.get("aa_name")

        if not isinstance(aa_id, int):
            logger.print("[ERROR] Invalid new_residue aa_id in clean report.")
            return None
        if not isinstance(aa_name, str) or aa_name.strip() == "":
            logger.print("[ERROR] Invalid new_residue aa_name in clean report.")
            return None

        new_residue_list.append({
            "aa_id": aa_id,
            "aa_name": aa_name.strip().upper(),
        })

    new_residue_list.sort(key=lambda x: x["aa_id"])
    return new_residue_list


def get_supported_output_type(data: Dict[str, Any], logger: Logger) -> str | None:
    output_type = data.get("output_type")
    if not isinstance(output_type, str):
        logger.print("[ERROR] Missing or invalid output_type.")
        return None
    if output_type not in SUPPORTED_OUTPUT_TYPES:
        logger.print(f"[ERROR] Unsupported output_type for integrate: {output_type}")
        return None
    return output_type


def validate_report_by_type(data: Dict[str, Any], logger: Logger) -> bool:
    output_type = get_supported_output_type(data, logger)
    if output_type is None:
        return False

    if output_type == "enzywizard_clean":
        return validate_clean_report(data, logger)
    if output_type == "enzywizard_aaprops":
        return validate_aaprops_report(data, logger)
    if output_type == "enzywizard_hydrocluster":
        return validate_hydrocluster_report(data, logger)
    if output_type == "enzywizard_energy":
        return validate_energy_report(data, logger)
    if output_type == "enzywizard_flexibility":
        return validate_flexibility_report(data, logger)
    if output_type == "enzywizard_disorder":
        return validate_disorder_report(data, logger)
    if output_type == "enzywizard_conservation":
        return validate_conservation_report(data, logger)
    if output_type == "enzywizard_embedding":
        return validate_embedding_report(data, logger)
    if output_type == "enzywizard_pocket":
        return validate_pocket_report(data, logger)
    if output_type == "enzywizard_substrate":
        return validate_substrate_report(data, logger)
    if output_type == "enzywizard_dock":
        return validate_dock_report(data, logger)
    if output_type == "enzywizard_interaction":
        return validate_interaction_report(data, logger)

    logger.print("[ERROR] Unsupported report type.")
    return False


def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _is_number_list(x: Any, min_len: int | None = None) -> bool:
    if not isinstance(x, list):
        return False
    if min_len is not None and len(x) < min_len:
        return False
    return all(_is_number(v) for v in x)


def _is_int01_list(x: Any) -> bool:
    if not isinstance(x, list):
        return False
    return all(isinstance(v, int) and v in {0, 1} for v in x)


def _validate_residue_id_name(item: Dict[str, Any], logger: Logger) -> bool:
    if not isinstance(item.get("aa_id"), int):
        logger.print("[ERROR] Invalid aa_id.")
        return False
    if not isinstance(item.get("aa_name"), str) or item["aa_name"].strip() == "":
        logger.print("[ERROR] Invalid aa_name.")
        return False
    return True


def validate_clean_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_clean":
        logger.print("[ERROR] clean report output_type mismatch.")
        return False

    mapping = data.get("amino_acid_mapping_old_to_new")
    stats = data.get("clean_statistics")

    if not isinstance(mapping, list):
        logger.print("[ERROR] Invalid amino_acid_mapping_old_to_new.")
        return False
    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid clean_statistics.")
        return False

    required_stat_keys = [
        "changed_resname",
        "removed_nonstd",
        "removed_missing_bb",
        "removed_missing_heavy_atoms",
        "removed_unexpected_heavy_atoms",
        "removed_bad_occ",
        "removed_inscodes",
        "kept_residues",
    ]
    for k in required_stat_keys:
        if not isinstance(stats.get(k), int):
            logger.print(f"[ERROR] Invalid clean_statistics field: {k}")
            return False

    for item in mapping:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid residue mapping entry.")
            return False
        old_residue = item.get("old_residue")
        new_residue = item.get("new_residue")
        if not isinstance(old_residue, dict) or not isinstance(new_residue, dict):
            logger.print("[ERROR] Invalid residue mapping structure.")
            return False
        for residue_item in [old_residue, new_residue]:
            if not _validate_residue_id_name(residue_item, logger):
                return False
            if not isinstance(residue_item.get("hydrogen_atom_count"), int):
                logger.print("[ERROR] Invalid hydrogen_atom_count.")
                return False

    return True


def validate_aaprops_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_aaprops":
        logger.print("[ERROR] aaprops report output_type mismatch.")
        return False

    aa_props = data.get("aa_props")
    stats = data.get("aa_props_statistics")

    if not isinstance(aa_props, list):
        logger.print("[ERROR] Invalid aa_props.")
        return False
    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid aa_props_statistics.")
        return False

    for key in ["aa_name_statistics", "aa_class_statistics", "aa_ss_statistics"]:
        if not isinstance(stats.get(key), dict):
            logger.print(f"[ERROR] Invalid aa_props_statistics field: {key}")
            return False

    for item in aa_props:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid aa_props entry.")
            return False

        required_int_fields = ["aa_id"]
        required_str_fields = ["aa_name", "aa_class", "aa_ss"]
        required_list_int_fields = ["aa_name_one_hot", "aa_class_one_hot", "aa_ss_one_hot"]
        required_float_fields = [
            "aa_rsa", "aa_phi", "aa_psi", "aa_net_charge", "aa_pka", "aa_volume",
            "aa_hydrophobicity", "aa_molecular_weight", "aa_pi"
        ]

        for k in required_int_fields:
            if not isinstance(item.get(k), int):
                logger.print(f"[ERROR] Invalid aaprops field: {k}")
                return False
        for k in required_str_fields:
            if not isinstance(item.get(k), str):
                logger.print(f"[ERROR] Invalid aaprops field: {k}")
                return False
        for k in required_list_int_fields:
            if not _is_int01_list(item.get(k)):
                logger.print(f"[ERROR] Invalid one-hot field: {k}")
                return False
        for k in required_float_fields:
            if not _is_number(item.get(k)):
                logger.print(f"[ERROR] Invalid numeric field: {k}")
                return False
        if not _is_number_list(item.get("aa_coord"), min_len=3):
            logger.print("[ERROR] Invalid aa_coord.")
            return False

    return True


def validate_hydrocluster_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_hydrocluster":
        logger.print("[ERROR] hydrocluster report output_type mismatch.")
        return False

    stats = data.get("hydrophobic_cluster_statistics")
    clusters = data.get("hydrophobic_cluster")

    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid hydrophobic_cluster_statistics.")
        return False
    if not isinstance(stats.get("cluster_num"), int):
        logger.print("[ERROR] Invalid hydrophobic_cluster_statistics.cluster_num.")
        return False
    if not _is_number(stats.get("max_cluster_area")):
        logger.print("[ERROR] Invalid hydrophobic_cluster_statistics.max_cluster_area.")
        return False
    if not _is_number(stats.get("total_cluster_area")):
        logger.print("[ERROR] Invalid hydrophobic_cluster_statistics.total_cluster_area.")
        return False

    if not isinstance(clusters, list):
        logger.print("[ERROR] Invalid hydrophobic_cluster.")
        return False

    for item in clusters:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid hydrophobic_cluster entry.")
            return False
        if not _is_number(item.get("area")):
            logger.print("[ERROR] Invalid hydrophobic cluster area.")
            return False
        residues = item.get("residues")
        if not isinstance(residues, list):
            logger.print("[ERROR] Invalid hydrophobic cluster residues.")
            return False
        for residue in residues:
            if not isinstance(residue, dict) or not _validate_residue_id_name(residue, logger):
                return False

    return True


def validate_energy_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_energy":
        logger.print("[ERROR] energy report output_type mismatch.")
        return False

    terms = data.get("energy_terms")
    if not isinstance(terms, dict):
        logger.print("[ERROR] Invalid energy_terms.")
        return False

    keys = [
        "total_potential_energy",
        "harmonic_bond_force",
        "harmonic_angle_force",
        "custom_bond_force",
        "custom_torsion_force",
        "custom_nonbonded_force",
        "nonbonded_force",
        "periodic_torsion_force",
        "cmap_torsion_force",
    ]
    for k in keys:
        if not _is_number(terms.get(k)):
            logger.print(f"[ERROR] Invalid energy term: {k}")
            return False

    return True


def validate_flexibility_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_flexibility":
        logger.print("[ERROR] flexibility report output_type mismatch.")
        return False

    entries = data.get("protein_rmsf")
    if not isinstance(entries, list):
        logger.print("[ERROR] Invalid protein_rmsf.")
        return False

    for item in entries:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid protein_rmsf entry.")
            return False
        if not _validate_residue_id_name(item, logger):
            return False
        if not _is_number(item.get("rmsf")):
            logger.print("[ERROR] Invalid rmsf.")
            return False

    return True


def validate_disorder_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_disorder":
        logger.print("[ERROR] disorder report output_type mismatch.")
        return False

    stats = data.get("disorder_region_statistics")
    regions = data.get("disorder_regions")

    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid disorder_region_statistics.")
        return False
    if not isinstance(stats.get("region_num"), int):
        logger.print("[ERROR] Invalid disorder_region_statistics.region_num.")
        return False
    if not isinstance(stats.get("max_region_length"), int):
        logger.print("[ERROR] Invalid disorder_region_statistics.max_region_length.")
        return False
    if not isinstance(stats.get("total_region_length"), int):
        logger.print("[ERROR] Invalid disorder_region_statistics.total_region_length.")
        return False

    if not isinstance(regions, list):
        logger.print("[ERROR] Invalid disorder_regions.")
        return False

    for item in regions:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid disorder region entry.")
            return False
        if not isinstance(item.get("length"), int):
            logger.print("[ERROR] Invalid disorder region length.")
            return False
        residues = item.get("residues")
        if not isinstance(residues, list):
            logger.print("[ERROR] Invalid disorder region residues.")
            return False
        for residue in residues:
            if not isinstance(residue, dict) or not _validate_residue_id_name(residue, logger):
                return False

    return True


def validate_conservation_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_conservation":
        logger.print("[ERROR] conservation report output_type mismatch.")
        return False

    entries = data.get("conservation_scores")
    if not isinstance(entries, list):
        logger.print("[ERROR] Invalid conservation_scores.")
        return False

    for item in entries:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid conservation entry.")
            return False
        if not _validate_residue_id_name(item, logger):
            return False
        for key in ["hmm_emission_log_score", "emission_probability", "conservation_score"]:
            if not _is_number(item.get(key)):
                logger.print(f"[ERROR] Invalid conservation field: {key}")
                return False

    return True


def validate_embedding_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_embedding":
        logger.print("[ERROR] embedding report output_type mismatch.")
        return False

    entries = data.get("embeddings")
    if not isinstance(entries, list):
        logger.print("[ERROR] Invalid embeddings.")
        return False

    for item in entries:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid embedding entry.")
            return False
        if not _validate_residue_id_name(item, logger):
            return False
        if not _is_number_list(item.get("embedding"), min_len=1):
            logger.print("[ERROR] Invalid embedding vector.")
            return False

    return True


def validate_pocket_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_pocket":
        logger.print("[ERROR] pocket report output_type mismatch.")
        return False

    stats = data.get("pocket_region_statistics")
    regions = data.get("pocket_regions")

    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid pocket_region_statistics.")
        return False
    if not isinstance(stats.get("pocket_num"), int):
        logger.print("[ERROR] Invalid pocket_region_statistics.pocket_num.")
        return False
    if not _is_number(stats.get("max_pocket_volume")):
        logger.print("[ERROR] Invalid pocket_region_statistics.max_pocket_volume.")
        return False
    if not _is_number(stats.get("total_pocket_volume")):
        logger.print("[ERROR] Invalid pocket_region_statistics.total_pocket_volume.")
        return False

    if not isinstance(regions, list):
        logger.print("[ERROR] Invalid pocket_regions.")
        return False

    for item in regions:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid pocket region entry.")
            return False
        if not _is_number(item.get("volume")):
            logger.print("[ERROR] Invalid pocket volume.")
            return False
        if not isinstance(item.get("n_spheres"), int):
            logger.print("[ERROR] Invalid pocket n_spheres.")
            return False
        if not _is_number_list(item.get("pocket_center_coord"), min_len=3):
            logger.print("[ERROR] Invalid pocket_center_coord.")
            return False
        if not _is_number_list(item.get("pocket_box_boundaries"), min_len=3):
            logger.print("[ERROR] Invalid pocket_box_boundaries.")
            return False
        residues = item.get("residues")
        if not isinstance(residues, list):
            logger.print("[ERROR] Invalid pocket residues.")
            return False
        for residue in residues:
            if not isinstance(residue, dict) or not _validate_residue_id_name(residue, logger):
                return False

    return True


def validate_substrate_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_substrate":
        logger.print("[ERROR] substrate report output_type mismatch.")
        return False

    substrates = data.get("substrates")
    if not isinstance(substrates, list):
        logger.print("[ERROR] Invalid substrates.")
        return False

    for item in substrates:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid substrate entry.")
            return False
        for key in ["substrate_name", "smiles"]:
            if not isinstance(item.get(key), str):
                logger.print(f"[ERROR] Invalid substrate field: {key}")
                return False
        if not _is_int01_list(item.get("fingerprint")):
            logger.print("[ERROR] Invalid fingerprint.")
            return False
        if not isinstance(item.get("num_atoms"), int):
            logger.print("[ERROR] Invalid num_atoms.")
            return False
        for key in ["mol_weight", "logp"]:
            if not _is_number(item.get(key)):
                logger.print(f"[ERROR] Invalid substrate field: {key}")
                return False

        structures = item.get("structures")
        if not isinstance(structures, list):
            logger.print("[ERROR] Invalid substrate structures.")
            return False
        for s in structures:
            if not isinstance(s, dict):
                logger.print("[ERROR] Invalid substrate structure entry.")
                return False
            if not isinstance(s.get("structure_name"), str):
                logger.print("[ERROR] Invalid structure_name.")
                return False
            if not _is_number(s.get("structure_energy")):
                logger.print("[ERROR] Invalid structure_energy.")
                return False

    return True


def validate_dock_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_dock":
        logger.print("[ERROR] dock report output_type mismatch.")
        return False

    result = data.get("docked_result")
    if not isinstance(result, dict):
        logger.print("[ERROR] Invalid docked_result.")
        return False

    for key in ["complex_name", "substrate_names"]:
        if not isinstance(result.get(key), str):
            logger.print(f"[ERROR] Invalid dock field: {key}")
            return False

    if not _is_number(result.get("docking_score")):
        logger.print("[ERROR] Invalid docking_score.")
        return False
    if not _is_number_list(result.get("docking_box_center"), min_len=3):
        logger.print("[ERROR] Invalid docking_box_center.")
        return False
    if not _is_number_list(result.get("docking_box_size"), min_len=3):
        logger.print("[ERROR] Invalid docking_box_size.")
        return False

    docked_substrates = result.get("docked_substrates")
    if not isinstance(docked_substrates, list):
        logger.print("[ERROR] Invalid docked_substrates.")
        return False

    for item in docked_substrates:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid docked substrate entry.")
            return False
        if not isinstance(item.get("substrate_name"), str):
            logger.print("[ERROR] Invalid docked substrate_name.")
            return False
        if not isinstance(item.get("conformation_name"), str):
            logger.print("[ERROR] Invalid conformation_name.")
            return False
        if not _is_number_list(item.get("docked_center_coord"), min_len=3):
            logger.print("[ERROR] Invalid docked_center_coord.")
            return False

    return True


def validate_interaction_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("output_type") != "enzywizard_interaction":
        logger.print("[ERROR] interaction report output_type mismatch.")
        return False

    interactions = data.get("interactions")
    stats = data.get("interactions_statistics")

    if not isinstance(interactions, list):
        logger.print("[ERROR] Invalid interactions.")
        return False
    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid interactions_statistics.")
        return False

    for block_name in ["overall", "intra_protein", "protein_substrate"]:
        block = stats.get(block_name)
        if not isinstance(block, dict):
            logger.print(f"[ERROR] Invalid interaction statistics block: {block_name}")
            return False
        for sub_key in ["count", "unique_pair_count"]:
            sub_block = block.get(sub_key)
            if not isinstance(sub_block, dict):
                logger.print(f"[ERROR] Invalid interaction statistics sub-block: {block_name}.{sub_key}")
                return False
            for interaction_name in INTERACTION_ORDER:
                if not isinstance(sub_block.get(interaction_name), int):
                    logger.print(f"[ERROR] Invalid interaction statistic: {block_name}.{sub_key}.{interaction_name}")
                    return False

    for item in interactions:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid interaction entry.")
            return False

        interaction = item.get("interaction")
        node1 = item.get("node1")
        node2 = item.get("node2")

        if interaction not in INTERACTION_ORDER:
            logger.print("[ERROR] Invalid interaction type.")
            return False
        if not isinstance(node1, dict) or not isinstance(node2, dict):
            logger.print("[ERROR] Invalid interaction node.")
            return False
        if not validate_interaction_node(node1, logger):
            return False
        if not validate_interaction_node(node2, logger):
            return False

    return True


def validate_interaction_node(node: Dict[str, Any], logger: Logger) -> bool:
    node_type = node.get("node_type")
    if node_type == "amino_acid":
        if not isinstance(node.get("aa_index"), int):
            logger.print("[ERROR] Invalid amino_acid node aa_index.")
            return False
        if not isinstance(node.get("aa_name"), str):
            logger.print("[ERROR] Invalid amino_acid node aa_name.")
            return False
        return True

    if node_type == "substrate":
        if not isinstance(node.get("substrate_index"), int):
            logger.print("[ERROR] Invalid substrate node substrate_index.")
            return False
        if not isinstance(node.get("substrate_name"), str):
            logger.print("[ERROR] Invalid substrate node substrate_name.")
            return False
        return True

    logger.print("[ERROR] Invalid interaction node_type.")
    return False