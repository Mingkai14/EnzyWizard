from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List
import re

from rdkit import Chem

from ..utils.logging_utils import Logger


def load_sdf_mol_3d(sdf_path: str | Path, logger: Logger) -> Chem.Mol | None:
    try:
        sdf_path = Path(sdf_path)

        if not sdf_path.exists() or sdf_path.stat().st_size <= 0:
            logger.print("[ERROR] Invalid input SDF file.")
            return None

        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        if supplier is None or len(supplier) == 0:
            logger.print("[ERROR] Failed to load SDF file.")
            return None

        mol = supplier[0]
        if mol is None:
            logger.print("[ERROR] Failed to parse Mol from SDF file.")
            return None

        if mol.GetNumConformers() <= 0:
            logger.print("[ERROR] Input SDF does not contain 3D coordinates.")
            return None

        return mol

    except Exception:
        logger.print("[ERROR] Failed to read Mol(3D) from SDF file.")
        return None


def get_sdf_atom_info(
    sdf_path: str | Path,
    logger: Logger,
) -> Dict[str, Any] | None:
    mol = load_sdf_mol_3d(sdf_path, logger)
    if mol is None:
        return None

    try:
        atom_info_list: List[Dict[str, Any]] = []

        for atom in mol.GetAtoms():
            atom_info_list.append(
                {
                    "atom_index": int(atom.GetIdx() + 1),
                    "atom_name": str(atom.GetSymbol()).upper(),
                }
            )

        return {
            "atom_count": int(mol.GetNumAtoms()),
            "atom_info_list": atom_info_list,
        }

    except Exception:
        logger.print("[ERROR] Failed to extract atom information from SDF file.")
        return None


def get_pdbqt_atom_info_from_lines(lines: List[str]) -> List[Dict[str, Any]] | None:
    try:
        atom_info_list: List[Dict[str, Any]] = []

        for line in lines:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            if len(line) < 54:
                return None

            try:
                pdbqt_atom_index = int(line[6:11].strip())
                pdbqt_atom_name = line[12:16].strip().upper()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except Exception:
                return None

            atom_info_list.append(
                {
                    "pdbqt_atom_index": pdbqt_atom_index,
                    "pdbqt_atom_name": pdbqt_atom_name,
                    "x": x,
                    "y": y,
                    "z": z,
                }
            )

        return atom_info_list

    except Exception:
        return None


def get_pdbqt_index_mapping(
    pdbqt_path: str | Path,
    logger: Logger,
) -> List[Dict[str, Any]] | None:
    try:
        pdbqt_path = Path(pdbqt_path)

        if not pdbqt_path.exists() or pdbqt_path.stat().st_size <= 0:
            logger.print("[ERROR] Invalid input PDBQT file.")
            return None

        lines = pdbqt_path.read_text(encoding="utf-8", errors="replace").splitlines()

        atom_info_list = get_pdbqt_atom_info_from_lines(lines)
        if atom_info_list is None or len(atom_info_list) == 0:
            logger.print("[ERROR] Failed to read atom information from PDBQT file.")
            return None

        pdbqt_index_to_name: Dict[int, str] = {}
        for item in atom_info_list:
            pdbqt_index_to_name[int(item["pdbqt_atom_index"])] = str(item["pdbqt_atom_name"])

        mapping_numbers: List[int] = []
        for line in lines:
            upper_line = line.upper()
            if "REMARK" not in upper_line:
                continue
            if not upper_line.startswith("REMARK"):
                continue
            if "INDEX MAP" not in upper_line:
                continue

            nums = re.findall(r"\d+", line)
            mapping_numbers.extend(int(x) for x in nums)

        if len(mapping_numbers) < 2 or len(mapping_numbers) % 2 != 0:
            logger.print("[ERROR] Invalid INDEX MAP format.")
            return None


        mapping_info_list: List[Dict[str, Any]] = []

        for i in range(0, len(mapping_numbers), 2):
            original_atom_index = int(mapping_numbers[i])
            pdbqt_atom_index = int(mapping_numbers[i + 1])

            if pdbqt_atom_index not in pdbqt_index_to_name:
                logger.print("[ERROR] PDBQT atom index in mapping not found in atom records.")
                return None

            mapping_info_list.append(
                {
                    "original_atom_index": original_atom_index,
                    "original_atom_name": "",
                    "pdbqt_atom_index": pdbqt_atom_index,
                    "pdbqt_atom_name": pdbqt_index_to_name[pdbqt_atom_index],
                }
            )

        mapping_info_list.sort(key=lambda x: int(x["pdbqt_atom_index"]))

        return mapping_info_list

    except Exception:
        logger.print("[ERROR] Failed to parse index mapping from PDBQT file.")
        return None


def get_pose_ligand_block_list(
    pose_string: str,
    logger: Logger,
) -> List[List[str]] | None:
    if not isinstance(pose_string, str) or len(pose_string.strip()) == 0:
        logger.print("[ERROR] Invalid pose string.")
        return None

    try:
        lines = pose_string.splitlines()
        ligand_block_list: List[List[str]] = []
        current_block: List[str] = []
        in_ligand_block = False

        for line in lines:
            stripped = line.strip()
            upper_line = stripped.upper()

            if upper_line == "ROOT":
                if in_ligand_block:
                    logger.print("[ERROR] Nested ROOT found.")
                    return None
                current_block = [line]
                in_ligand_block = True
                continue

            if in_ligand_block:
                current_block.append(line)

                if upper_line.startswith("TORSDOF"):
                    ligand_block_list.append(current_block)
                    current_block = []
                    in_ligand_block = False

        if in_ligand_block:
            logger.print("[ERROR] Incomplete ligand block (missing TORSDOF).")
            return None

        if len(ligand_block_list) == 0:
            logger.print("[ERROR] No ligand blocks found in pose string.")
            return None

        return ligand_block_list

    except Exception:
        logger.print("[ERROR] Failed to split pose string into ligand blocks.")
        return None


def get_pose_for_substrate_atom_info(
    substrate_name: str,
    ligand_order_index: int,
    pose_string: str,
    original_atom_count: int,
    original_atom_info_list: List[Dict[str, Any]],
    mapping_info_list: List[Dict[str, Any]],
    logger: Logger,
) -> Dict[str, Any] | None:
    if not substrate_name:
        logger.print("[ERROR] substrate_name is empty.")
        return None

    if ligand_order_index < 0:
        logger.print("[ERROR] ligand_order_index must be non-negative.")
        return None

    if original_atom_count <= 0:
        logger.print("[ERROR] original_atom_count must be positive.")
        return None

    if not isinstance(original_atom_info_list, list) or len(original_atom_info_list) != original_atom_count:
        logger.print("[ERROR] Invalid original_atom_info_list.")
        return None

    if not isinstance(mapping_info_list, list) or len(mapping_info_list) != original_atom_count:
        logger.print("[ERROR] Invalid mapping_info_list.")
        return None

    ligand_block_list = get_pose_ligand_block_list(pose_string, logger)
    if ligand_block_list is None:
        return None

    if ligand_order_index >= len(ligand_block_list):
        logger.print(f"[ERROR] ligand_order_index {ligand_order_index} is out of range.")
        return None

    try:
        original_index_to_name: Dict[int, str] = {}
        for item in original_atom_info_list:
            atom_index = int(item.get("atom_index", 0))
            atom_name = str(item.get("atom_name", "")).upper()

            if atom_index <= 0 or not atom_name:
                logger.print("[ERROR] Invalid original atom information.")
                return None

            original_index_to_name[atom_index] = atom_name

        enriched_mapping_info_list: List[Dict[str, Any]] = []
        for item in mapping_info_list:
            original_atom_index = int(item.get("original_atom_index", 0))
            pdbqt_atom_index = int(item.get("pdbqt_atom_index", 0))
            pdbqt_atom_name = str(item.get("pdbqt_atom_name", "")).upper()

            if (
                original_atom_index <= 0
                or pdbqt_atom_index <= 0
                or not pdbqt_atom_name
                or original_atom_index not in original_index_to_name
            ):
                logger.print("[ERROR] Invalid mapping information.")
                return None

            enriched_mapping_info_list.append(
                {
                    "original_atom_index": original_atom_index,
                    "original_atom_name": original_index_to_name[original_atom_index],
                    "pdbqt_atom_index": pdbqt_atom_index,
                    "pdbqt_atom_name": pdbqt_atom_name,
                }
            )

        if len(enriched_mapping_info_list) != original_atom_count:
            logger.print("[ERROR] Mapping size mismatch with original atoms.")
            return None

        enriched_mapping_info_list.sort(key=lambda x: int(x["pdbqt_atom_index"]))

        expected_index_set = set(
            int(item["pdbqt_atom_index"]) for item in enriched_mapping_info_list
        )

        expected_pdbqt_atom_name_list = [
            str(item["pdbqt_atom_name"]).upper()
            for item in enriched_mapping_info_list
        ]

        block_lines = ligand_block_list[ligand_order_index]
        matched_atom_info_list = get_pdbqt_atom_info_from_lines(block_lines)
        if matched_atom_info_list is None:
            logger.print("[ERROR] Failed to parse atom information from pose ligand block.")
            return None

        if len(matched_atom_info_list) != original_atom_count:
            logger.print(f"[ERROR] Atom count mismatch for substrate: {substrate_name}")
            return None

        matched_atom_info_list.sort(key=lambda x: int(x["pdbqt_atom_index"]))

        pose_index_set = set(
            int(item["pdbqt_atom_index"]) for item in matched_atom_info_list
        )
        if pose_index_set != expected_index_set:
            logger.print(f"[ERROR] PDBQT atom index mismatch for substrate: {substrate_name}")
            return None

        pose_pdbqt_atom_name_list = [
            str(item["pdbqt_atom_name"]).upper()
            for item in matched_atom_info_list
        ]
        if pose_pdbqt_atom_name_list != expected_pdbqt_atom_name_list:
            logger.print(f"[ERROR] PDBQT atom name mismatch for substrate: {substrate_name}")
            return None

        pdbqt_index_to_pose_atom: Dict[int, Dict[str, Any]] = {}
        for item in matched_atom_info_list:
            pdbqt_index_to_pose_atom[int(item["pdbqt_atom_index"])] = item

        docked_atom_info_list: List[Dict[str, Any]] = []
        for item in enriched_mapping_info_list:
            original_atom_index = int(item["original_atom_index"])
            pdbqt_atom_index = int(item["pdbqt_atom_index"])

            if pdbqt_atom_index not in pdbqt_index_to_pose_atom:
                logger.print("[ERROR] Matched ligand block is inconsistent with mapping.")
                return None

            pose_atom = pdbqt_index_to_pose_atom[pdbqt_atom_index]

            docked_atom_info_list.append(
                {
                    "original_atom_index": original_atom_index,
                    "original_atom_name": str(item["original_atom_name"]).upper(),
                    "pdbqt_atom_index": pdbqt_atom_index,
                    "pdbqt_atom_name": str(item["pdbqt_atom_name"]).upper(),
                    "x": float(pose_atom["x"]),
                    "y": float(pose_atom["y"]),
                    "z": float(pose_atom["z"]),
                }
            )

        docked_atom_info_list.sort(key=lambda x: int(x["original_atom_index"]))

        return {
            "substrate_name": substrate_name,
            "atom_info_list": docked_atom_info_list,
        }

    except Exception:
        logger.print(f"[ERROR] Failed to parse pose for substrate: {substrate_name}")
        return None


def write_docked_sdf_from_atom_info(
    original_mol_3d: Chem.Mol,
    docked_atom_info_list: List[Dict[str, Any]],
    sdf_path: str | Path,
    logger: Logger,
) -> bool:
    if original_mol_3d is None or original_mol_3d.GetNumConformers() <= 0:
        logger.print("[ERROR] Invalid original Mol(3D).")
        return False

    if not isinstance(docked_atom_info_list, list) or len(docked_atom_info_list) != original_mol_3d.GetNumAtoms():
        logger.print("[ERROR] Invalid docked_atom_info_list.")
        return False

    try:
        sdf_path = Path(sdf_path)
        sdf_path.parent.mkdir(parents=True, exist_ok=True)

        mol = Chem.Mol(original_mol_3d)
        conf = mol.GetConformer()

        for item in docked_atom_info_list:
            original_atom_index = int(item.get("original_atom_index", 0))
            x = float(item.get("x", 0.0))
            y = float(item.get("y", 0.0))
            z = float(item.get("z", 0.0))

            if original_atom_index <= 0 or original_atom_index > mol.GetNumAtoms():
                logger.print("[ERROR] Invalid original atom index in docked_atom_info_list.")
                return False

            atom = mol.GetAtomWithIdx(original_atom_index - 1)
            original_atom_name = str(item.get("original_atom_name", "")).upper()
            if original_atom_name and atom.GetSymbol().upper() != original_atom_name:
                logger.print("[ERROR] Atom name mismatch when writing docked SDF.")
                return False

            conf.SetAtomPosition(original_atom_index - 1, (x, y, z))

        writer = Chem.SDWriter(str(sdf_path))
        conf_id = conf.GetId()
        writer.write(mol, confId=conf_id)
        writer.close()

        if not sdf_path.exists() or sdf_path.stat().st_size <= 0:
            logger.print("[ERROR] Failed to save docked SDF file.")
            return False

        return True

    except Exception:
        logger.print("[ERROR] Failed to write docked atom information to SDF file.")
        return False