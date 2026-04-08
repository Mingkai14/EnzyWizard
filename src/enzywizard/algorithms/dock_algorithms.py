from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List
import tempfile
from itertools import product

from vina import Vina
from ..utils.dock_utils import split_vina_pose_string, get_sdf_atom_info_from_mol, get_pdbqt_index_mapping, get_pose_for_substrate_atom_info
from ..utils.IO_utils import load_sdf_mol_3d, write_protein_pdbqt,write_substrate_pdbqt_from_sdf, write_docked_complex_cif, write_docked_sdf_from_atom_info
from Bio.PDB.Structure import Structure
from ..utils.logging_utils import Logger
from ..algorithms.pocket_algorithms import compute_pockets
from ..utils.structure_utils import get_structure_box
from ..utils.dock_utils import get_substrate_sdf_path_group_dict, compute_ligand_centroid


def dock_multiple_ligands_with_vina(
    protein_pdbqt_path: str | Path,
    ligand_pdbqt_path_list: List[str | Path],
    ligand_sdf_path_list: List[str | Path],
    ligand_name_list: List[str],
    box_center_list: List[float],
    box_size_list: List[float],
    logger: Logger,
    exhaustiveness: int = 16,
    cpu: int = 0,
    max_pose_read_num: int = 3
) -> List[Dict[str, Any]] | None:

    if not isinstance(protein_pdbqt_path, (str, Path)):
        logger.print("[ERROR] protein_pdbqt_path must be a str or Path.")
        return None

    if not isinstance(ligand_pdbqt_path_list, list):
        logger.print("[ERROR] ligand_pdbqt_path_list must be a list.")
        return None

    if not isinstance(ligand_sdf_path_list, list):
        logger.print("[ERROR] ligand_sdf_path_list must be a list.")
        return None

    if not isinstance(ligand_name_list, list):
        logger.print("[ERROR] ligand_name_list must be a list.")
        return None

    if not isinstance(box_center_list, list):
        logger.print("[ERROR] box_center_list must be a list.")
        return None

    if not isinstance(box_size_list, list):
        logger.print("[ERROR] box_size_list must be a list.")
        return None

    if len(ligand_pdbqt_path_list) == 0:
        logger.print("[ERROR] ligand_pdbqt_path_list must not be empty.")
        return None

    if len(ligand_pdbqt_path_list) != len(ligand_sdf_path_list) or len(ligand_pdbqt_path_list) != len(ligand_name_list):
        logger.print("[ERROR] ligand_pdbqt_path_list, ligand_sdf_path_list, and ligand_name_list must have the same length.")
        return None

    if len(box_center_list) != 3:
        logger.print("[ERROR] box_center_list must contain exactly 3 values.")
        return None

    if len(box_size_list) != 3:
        logger.print("[ERROR] box_size_list must contain exactly 3 values.")
        return None

    if not isinstance(exhaustiveness, int) or exhaustiveness <= 0:
        logger.print("[ERROR] exhaustiveness must be a positive integer.")
        return None

    if not isinstance(cpu, int) or cpu < 0:
        logger.print("[ERROR] cpu must be a non-negative integer.")
        return None

    if not isinstance(max_pose_read_num, int) or max_pose_read_num <= 0:
        logger.print("[ERROR] max_pose_read_num must be a positive integer.")
        return None

    protein_pdbqt_path = Path(protein_pdbqt_path)

    if not protein_pdbqt_path.exists():
        logger.print(f"[ERROR] Protein PDBQT file does not exist: {protein_pdbqt_path}")
        return None

    try:
        box_center = [float(box_center_list[0]), float(box_center_list[1]), float(box_center_list[2])]
        box_size = [float(box_size_list[0]), float(box_size_list[1]), float(box_size_list[2])]
    except Exception:
        logger.print("[ERROR] box_center_list and box_size_list must contain numeric values.")
        return None

    if box_size[0] <= 0.0 or box_size[1] <= 0.0 or box_size[2] <= 0.0:
        logger.print("[ERROR] box_size_list values must be positive.")
        return None

    substrate_name_list: List[str] = []
    original_atom_info_list_list: List[List[Dict[str, Any]]] = []
    mapping_info_list_list: List[List[Dict[str, Any]]] = []
    ligand_pdbqt_str_list: List[str] = []

    try:
        for i in range(len(ligand_pdbqt_path_list)):
            ligand_pdbqt_path = Path(ligand_pdbqt_path_list[i])
            ligand_sdf_path = Path(ligand_sdf_path_list[i])
            ligand_name = str(ligand_name_list[i])

            if not ligand_pdbqt_path.exists():
                logger.print(f"[ERROR] Ligand PDBQT file does not exist: {ligand_pdbqt_path}")
                return None

            if not ligand_sdf_path.exists():
                logger.print(f"[ERROR] Ligand SDF file does not exist: {ligand_sdf_path}")
                return None

            pdbqt_stem = ligand_pdbqt_path.stem
            sdf_stem = ligand_sdf_path.stem
            name_stem = Path(ligand_name).stem

            if pdbqt_stem != sdf_stem or pdbqt_stem != name_stem:
                logger.print(
                    f"[ERROR] Stem mismatch at ligand index {i}: "
                    f"pdbqt_stem={pdbqt_stem}, sdf_stem={sdf_stem}, name_stem={name_stem}"
                )
                return None

            mol = load_sdf_mol_3d(ligand_sdf_path, logger)
            if mol is None:
                logger.print(f"[ERROR] Failed to load 3D SDF Mol: {ligand_sdf_path}")
                return None

            original_atom_info_list = get_sdf_atom_info_from_mol(mol, logger)
            if original_atom_info_list is None or len(original_atom_info_list) == 0:
                logger.print(f"[ERROR] Failed to read original atom info from SDF: {ligand_sdf_path}")
                return None

            mapping_info_list = get_pdbqt_index_mapping(ligand_pdbqt_path, logger)
            if mapping_info_list is None or len(mapping_info_list) == 0:
                logger.print(f"[ERROR] Failed to read PDBQT index mapping: {ligand_pdbqt_path}")
                return None

            substrate_name_list.append(pdbqt_stem)
            original_atom_info_list_list.append(original_atom_info_list)
            mapping_info_list_list.append(mapping_info_list)
            ligand_pdbqt_str_list.append(str(ligand_pdbqt_path))

    except Exception:
        logger.print("[ERROR] Failed during ligand input preparation.")
        return None

    try:
        v = Vina(cpu=cpu, verbosity=0)
        v.set_receptor(rigid_pdbqt_filename=str(protein_pdbqt_path))

        v.set_ligand_from_file(ligand_pdbqt_str_list)

        v.compute_vina_maps(center=box_center, box_size=box_size)
        v.dock(exhaustiveness=exhaustiveness)

    except Exception:
        logger.print("[ERROR] Vina docking failed.")
        return None

    try:
        pose_string = v.poses()
    except Exception:
        logger.print("[ERROR] Failed to read Vina poses.")
        return None


    try:
        energies = v.energies()
    except Exception:
        logger.print("[ERROR] Failed to read Vina energies.")
        return None

    if not pose_string and not energies:
        logger.print(f"[WARNING] Vina docking output is empty for: {';'.join(substrate_name_list)} in: {box_center} {box_size}.")
        return []

    if pose_string is None or energies is None:
        logger.print("[ERROR] Error during Vina outputting.")
        return None

    pose_string_list = split_vina_pose_string(pose_string,logger)

    if len(pose_string_list) == 0:
        logger.print("[ERROR] No docking pose was parsed from Vina output.")
        return None

    try:
        energy_value_list: List[float] = []

        for item in energies:
            try:
                energy_value_list.append(float(item[0]))
            except Exception:
                energy_value_list.append(float(item))

    except Exception:
        logger.print("[ERROR] Failed to parse Vina energies.")
        return None

    if len(pose_string_list) != len(energy_value_list):
        logger.print("[ERROR] Vina output poses are not equal with energies.")
        return None

    pose_num = min(max_pose_read_num, len(pose_string_list), len(energy_value_list))

    if pose_num <= 0:
        logger.print("[ERROR] No valid docking poses available after pose/energy matching.")
        return None

    result_list: List[Dict[str, Any]] = []
    joined_substrate_names = ";".join(substrate_name_list)

    try:
        for pose_index in range(pose_num):
            current_pose_string = pose_string_list[pose_index]
            current_energy = float(energy_value_list[pose_index])

            docked_substrate_info_list: List[Dict[str, Any]] = []

            for ligand_order_index in range(len(substrate_name_list)):
                substrate_name = substrate_name_list[ligand_order_index]
                original_atom_info_list = original_atom_info_list_list[ligand_order_index]
                mapping_info_list = mapping_info_list_list[ligand_order_index]

                parsed_substrate_info = get_pose_for_substrate_atom_info(
                    substrate_name=substrate_name,
                    ligand_order_index=ligand_order_index,
                    pose_string=current_pose_string,
                    original_atom_info_list=original_atom_info_list,
                    mapping_info_list=mapping_info_list,
                    logger=logger,
                )

                if parsed_substrate_info is None:
                    logger.print(
                        f"[ERROR] Failed to parse docked atom coordinates for substrate {substrate_name} in pose index {pose_index}."
                    )
                    return None

                parsed_substrate_info["source_sdf_path"] = str(Path(ligand_sdf_path_list[ligand_order_index]))
                docked_substrate_info_list.append(parsed_substrate_info)

            result_list.append(
                {
                    "substrate_names": joined_substrate_names,
                    "energy": current_energy,
                    "docked_substrate_info_list": docked_substrate_info_list,
                }
            )

        return result_list

    except Exception:
        logger.print("[ERROR] Failed to build docking result list.")
        return None

def dock_multiple_substrates_from_structure(
    struct: Structure,
    substrate_names: str,
    substrate_dir: str | Path,
    logger: Logger,
    max_docking_result_num: int = 3,
    max_docking_attempt_num: int = 20,
    max_pose_read_num: int = 1,
    exhaustiveness: int = 16,
    cpu: int = 0,
    min_rad: float = 1.8,
    max_rad: float = 6.2,
    min_volume: int = 50
) -> List[Dict[str, Any]] | None:
    if struct is None:
        logger.print("[ERROR] struct is None.")
        return None

    if not isinstance(substrate_names, str):
        logger.print("[ERROR] substrate_names must be a string.")
        return None

    if not substrate_names.strip():
        logger.print("[ERROR] substrate_names is empty.")
        return None

    if not isinstance(substrate_dir, (str, Path)):
        logger.print("[ERROR] substrate_dir must be a str or Path.")
        return None

    if not isinstance(max_docking_result_num, int) or max_docking_result_num <= 0:
        logger.print("[ERROR] max_docking_result_num must be a positive integer.")
        return None

    if not isinstance(max_docking_attempt_num, int) or max_docking_attempt_num <= 0:
        logger.print("[ERROR] max_docking_attempt_num must be a positive integer.")
        return None

    if not isinstance(max_pose_read_num, int) or max_pose_read_num <= 0:
        logger.print("[ERROR] max_pose_read_num must be a positive integer.")
        return None

    if not isinstance(exhaustiveness, int) or exhaustiveness <= 0:
        logger.print("[ERROR] exhaustiveness must be a positive integer.")
        return None

    if not isinstance(cpu, int) or cpu < 0:
        logger.print("[ERROR] cpu must be a non-negative integer.")
        return None

    if min_rad <= 0.0:
        logger.print("[ERROR] min_rad must be positive.")
        return None

    if max_rad <= 0.0:
        logger.print("[ERROR] max_rad must be positive.")
        return None

    if min_volume < 0:
        logger.print("[ERROR] min_volume must be non-negative.")
        return None

    if min_rad > max_rad:
        logger.print("[ERROR] min_rad must not be greater than max_rad.")
        return None

    substrate_dir = Path(substrate_dir)

    pocket_result_list = compute_pockets(
        struct=struct,
        logger=logger,
        min_rad=min_rad,
        max_rad=max_rad,
        min_volume=min_volume,
    )
    if pocket_result_list is None:
        logger.print("[ERROR] Failed to compute pockets.")
        return None

    structure_box_info = get_structure_box(struct=struct, logger=logger)
    if structure_box_info is None:
        logger.print("[ERROR] Failed to compute structure box.")
        return None

    grouped_result = get_substrate_sdf_path_group_dict(
        substrate_names=substrate_names,
        substrate_dir=substrate_dir,
        logger=logger,
    )
    if grouped_result is None:
        logger.print("[ERROR] Failed to group substrate SDF files.")
        return None

    substrate_name_list, substrate_sdf_path_group_dict = grouped_result

    if not isinstance(substrate_name_list, list) or len(substrate_name_list) == 0:
        logger.print("[ERROR] Invalid substrate_name_list.")
        return None

    if not isinstance(substrate_sdf_path_group_dict, dict) or len(substrate_sdf_path_group_dict) == 0:
        logger.print("[ERROR] Invalid substrate_sdf_path_group_dict.")
        return None

    substrate_sdf_path_list_list: List[List[str]] = []
    for substrate_name in substrate_name_list:
        sdf_path_list = substrate_sdf_path_group_dict.get(substrate_name, [])
        if not isinstance(sdf_path_list, list) or len(sdf_path_list) == 0:
            logger.print(f"[ERROR] No protomer SDF files found for substrate: {substrate_name}")
            return None
        substrate_sdf_path_list_list.append(sdf_path_list)

    box_info_list: List[Dict[str, Any]] = []

    for pocket_result in pocket_result_list:
        pocket_center_coord = pocket_result.get("pocket_center_coord", None)
        pocket_box_boundaries = pocket_result.get("pocket_box_boundaries", None)

        if (
            not isinstance(pocket_center_coord, list)
            or len(pocket_center_coord) != 3
            or not isinstance(pocket_box_boundaries, list)
            or len(pocket_box_boundaries) != 3
        ):
            logger.print("[ERROR] Invalid pocket box information.")
            return None

        box_info_list.append(
            {
                "box_center_list": [float(x) for x in pocket_center_coord],
                "box_size_list": [float(x) for x in pocket_box_boundaries],
            }
        )

    center_coord = structure_box_info.get("center_coord", None)
    box_boundaries = structure_box_info.get("box_boundaries", None)

    if (
        not isinstance(center_coord, list)
        or len(center_coord) != 3
        or not isinstance(box_boundaries, list)
        or len(box_boundaries) != 3
    ):
        logger.print("[ERROR] Invalid structure box information.")
        return None

    box_info_list.append(
        {
            "box_center_list": [float(x) for x in center_coord],
            "box_size_list": [float(x) for x in box_boundaries],
        }
    )

    result_list: List[Dict[str, Any]] = []
    docking_attempt_count = 0
    stop_docking = False

    try:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir_path = Path(tmp_dir)

            protein_pdbqt_path = tmp_dir_path / "protein.pdbqt"
            ok = write_protein_pdbqt(
                struct=struct,
                pdbqt_path=protein_pdbqt_path,
                logger=logger,
            )
            if not ok:
                logger.print("[ERROR] Failed to write protein PDBQT file.")
                return None

            sdf_to_pdbqt_path_dict: Dict[str, str] = {}

            for sdf_path_list in substrate_sdf_path_list_list:
                for sdf_path in sdf_path_list:
                    sdf_path_obj = Path(sdf_path)

                    if not sdf_path_obj.exists() or sdf_path_obj.stat().st_size <= 0:
                        logger.print(f"[ERROR] Invalid substrate SDF file: {sdf_path_obj}")
                        return None

                    pdbqt_path = tmp_dir_path / f"{sdf_path_obj.stem}.pdbqt"

                    ok = write_substrate_pdbqt_from_sdf(
                        sdf_path=sdf_path_obj,
                        pdbqt_path=pdbqt_path,
                        logger=logger,
                    )
                    if not ok:
                        logger.print(f"[ERROR] Failed to write substrate PDBQT file: {sdf_path_obj}")
                        return None

                    sdf_to_pdbqt_path_dict[str(sdf_path_obj)] = str(pdbqt_path)

            for selected_sdf_path_tuple in product(*substrate_sdf_path_list_list):
                if stop_docking:
                    break
                ligand_sdf_path_list: List[str] = [str(Path(x)) for x in selected_sdf_path_tuple]
                ligand_pdbqt_path_list: List[str] = []
                ligand_name_list: List[str] = []

                for sdf_path in ligand_sdf_path_list:
                    ligand_pdbqt_path = sdf_to_pdbqt_path_dict.get(str(Path(sdf_path)), "")
                    if not ligand_pdbqt_path:
                        logger.print(f"[ERROR] Missing cached ligand PDBQT path for: {sdf_path}")
                        return None

                    ligand_pdbqt_path_list.append(ligand_pdbqt_path)
                    ligand_name_list.append(Path(sdf_path).name)

                for box_info in box_info_list:
                    if stop_docking:
                        break
                    box_center_list = box_info["box_center_list"]
                    box_size_list = box_info["box_size_list"]

                    if docking_attempt_count >= max_docking_attempt_num:
                        logger.print("[WARNING] Maximum docking attempt count reached.")
                        stop_docking = True
                        break

                    docking_attempt_count += 1

                    docking_result_list = dock_multiple_ligands_with_vina(
                        protein_pdbqt_path=protein_pdbqt_path,
                        ligand_pdbqt_path_list=ligand_pdbqt_path_list,
                        ligand_sdf_path_list=ligand_sdf_path_list,
                        ligand_name_list=ligand_name_list,
                        box_center_list=box_center_list,
                        box_size_list=box_size_list,
                        exhaustiveness=exhaustiveness,
                        cpu=cpu,
                        max_pose_read_num=max_pose_read_num,
                        logger=logger,
                    )

                    if docking_result_list is None:
                        logger.print("[ERROR] dock_multiple_ligands_with_vina failed.")
                        return None

                    if len(docking_result_list) == 0:
                        continue

                    result_list.extend(docking_result_list)

                    if len(result_list) >= max_docking_result_num:
                        return result_list[:max_docking_result_num]

            if len(result_list) == 0:
                logger.print("[WARNING] No valid docking results were found for any substrate combination and docking box.")
                return []
            return result_list

    except Exception:
        logger.print("[ERROR] Failed to perform substrate docking workflow.")
        return None

def save_docking_results_and_generate_dock_report(
    docking_result_list: List[Dict[str, Any]],
    struct: Structure,
    protein_name: str,
    output_dir: str | Path,
    logger: Logger,
) -> Dict[str, Any] | None:

    if not isinstance(docking_result_list, list):
        logger.print("[ERROR] docking_result_list must be a list.")
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report_list: List[Dict[str, Any]] = []

    try:
        for i, docking_result in enumerate(docking_result_list):

            substrate_names = docking_result["substrate_names"]
            energy = float(docking_result["energy"])
            docked_list = docking_result["docked_substrate_info_list"]

            complex_name = f"docked_{protein_name}_{substrate_names}"

            # ===== 写 CIF =====
            cif_path = write_docked_complex_cif(
                docking_result=docking_result,
                struct=struct,
                protein_name=protein_name,
                output_dir=output_dir,
                logger=logger,
            )
            if cif_path is None:
                return None

            ligand_report_list: List[Dict[str, Any]] = []

            # ===== 写每个 ligand SDF =====
            for ligand in docked_list:

                substrate_name = ligand["substrate_name"]
                atom_info_list = ligand["atom_info_list"]
                source_sdf_path = ligand.get("source_sdf_path", "")

                mol = load_sdf_mol_3d(source_sdf_path, logger)
                if mol is None:
                    return None

                sdf_path = output_dir / f"docked_{substrate_name}.sdf"

                ok = write_docked_sdf_from_atom_info(
                    original_mol_3d=mol,
                    docked_atom_info_list=atom_info_list,
                    sdf_path=sdf_path,
                    logger=logger,
                )
                if not ok:
                    return None

                centroid = compute_ligand_centroid(atom_info_list, logger)
                if centroid is None:
                    return None

                ligand_report_list.append(
                    {
                        "substrate_name": substrate_name,
                        "centroid": centroid,
                    }
                )

            report_list.append(
                {
                    "complex_name": complex_name,
                    "substrate_names": substrate_names,
                    "energy": energy,
                    "ligands": ligand_report_list,
                }
            )

        return {
            "output_type": "enzywizard_dock",
            "docked_substrates": report_list,
        }

    except Exception:
        logger.print("[ERROR] Failed to save docking results.")
        return None