from __future__ import annotations
from Bio.PDB.Structure import Structure
from ..utils.logging_utils import Logger
from ..utils.structure_utils import get_single_chain, get_residues_by_chain
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from typing import Dict, List, Tuple, Optional, Union
from ..utils.clean_utils import (
    standardize_resname,
    choose_atom_altloc,
    clone_atom,
    normalize_atom_name,
    is_hydrogen_atom,
)
from ..resources.aa_resources import (
    AA3_STANDARD,
    AA3_REQUIRED_HEAVY_ATOMS,
    AA3_EXPECTED_HEAVY_ATOM_SET,
    AA3_ALLOWED_HEAVY_ATOM_SET_WITH_OXT,
    BACKBONE_REQUIRED_ATOMS,
)
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from openmm.app import Modeller, ForceField, PDBFile
import math
from ..utils.sequence_utils import normalize_aa_name_to_one_letter

def clean_structure_to_single_chain_A(struct: Structure, logger: Logger) -> Tuple[Structure, Dict[Tuple[int, str, str], Tuple[int, str, str]], Dict[str, int]] | None:
    old_chain = get_single_chain(struct, logger)
    if old_chain is None:
        return None

    new_struct = Structure(struct.id + "_cleaned")
    new_model = Model(0)
    new_chain = Chain("A")

    mapping_old_to_new: Dict[Tuple[int, str, str], Tuple[int, str, str]] = {}

    new_resseq = 0
    removed_nonstd = 0
    removed_missing_bb = 0
    removed_bad_occ = 0
    changed_resname = 0
    removed_inscodes = 0
    removed_missing_heavy_atoms = 0
    removed_unexpected_heavy_atoms = 0

    for res in old_chain.get_residues():
        hetflag, resseq, icode = res.id

        resname_old = res.get_resname().strip()
        resname_std = standardize_resname(resname_old)
        if resname_std != resname_old:
            changed_resname += 1

        if resname_std not in AA3_STANDARD:
            removed_nonstd += 1
            continue

        atoms_by_name: Dict[str, List[Atom]] = {}
        actual_heavy_atom_set = set()

        for atom in res.get_atoms():
            atom_name = normalize_atom_name(atom.get_name())
            atoms_by_name.setdefault(atom_name, []).append(atom)

            if not is_hydrogen_atom(atom):
                actual_heavy_atom_set.add(atom_name)

        if any(atom_name not in actual_heavy_atom_set for atom_name in BACKBONE_REQUIRED_ATOMS):
            removed_missing_bb += 1
            continue

        expected_heavy_atom_set = AA3_EXPECTED_HEAVY_ATOM_SET.get(resname_std)
        if expected_heavy_atom_set is None:
            removed_nonstd += 1
            continue

        if not expected_heavy_atom_set.issubset(actual_heavy_atom_set):
            removed_missing_heavy_atoms += 1
            continue

        allowed_heavy_atom_set = AA3_ALLOWED_HEAVY_ATOM_SET_WITH_OXT[resname_std]
        if not actual_heavy_atom_set.issubset(allowed_heavy_atom_set):
            removed_unexpected_heavy_atoms += 1
            continue

        bad_occ = False
        for atom_name in BACKBONE_REQUIRED_ATOMS:
            chosen = choose_atom_altloc(atoms_by_name[atom_name])
            occ = chosen.get_occupancy()
            if occ is not None and occ < 0:
                bad_occ = True
                break
        if bad_occ:
            removed_bad_occ += 1
            continue

        old_resseq = int(resseq)
        old_icode = str(icode)
        if old_icode.strip():
            removed_inscodes += 1

        new_resseq += 1
        new_icode = " "
        mapping_old_to_new[(old_resseq, resname_old, old_icode)] = (new_resseq, resname_std, new_icode)

        new_res = Residue((" ", new_resseq, " "), resname_std, res.get_segid())

        for atom_name, atom_list in atoms_by_name.items():
            chosen = choose_atom_altloc(atom_list)
            new_res.add(clone_atom(chosen))

        new_chain.add(new_res)

    new_model.add(new_chain)
    new_struct.add(new_model)

    stats = {
        "changed_resname": changed_resname,
        "removed_nonstd": removed_nonstd,
        "removed_missing_bb": removed_missing_bb,
        "removed_missing_heavy_atoms": removed_missing_heavy_atoms,
        "removed_unexpected_heavy_atoms": removed_unexpected_heavy_atoms,
        "removed_bad_occ": removed_bad_occ,
        "removed_inscodes": removed_inscodes,
        "kept_residues": new_resseq,
    }
    return new_struct, mapping_old_to_new, stats

def add_hydrogens_to_pdbfile(pdb: PDBFile, logger: Logger, pH:float = 7.0, force_field_file: str = "charmm36.xml") -> Modeller | None:

    try:
        ff = ForceField(force_field_file)

        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield=ff, pH=pH)

        return modeller

    except Exception as e:
        logger.print(f"[ERROR] Failed to add hydrogens: {e}")
        return None


def validate_clean_mapping_coordinates(structure: Structure,cleaned_structure: Structure,mapping_old_to_new: Dict[Tuple[int, str, str], Tuple[int, str, str]],logger: Logger) -> bool:
    tol: float = 1e-3

    old_chain = get_single_chain(structure, logger)
    if old_chain is None:
        return False

    new_chain = get_single_chain(cleaned_structure, logger)
    if new_chain is None:
        return False

    old_res_list = get_residues_by_chain(old_chain, logger)
    if old_res_list is None:
        return False

    new_res_list = get_residues_by_chain(new_chain, logger)
    if new_res_list is None:
        return False

    # 构建 dict
    old_dict = {}
    for (res_id, resname, coord) in old_res_list:
        _, resseq, icode = res_id
        old_dict[(int(resseq), resname, str(icode))] = coord

    new_dict = {}
    for (res_id, resname, coord) in new_res_list:
        _, resseq, icode = res_id
        new_dict[(int(resseq), resname, str(icode))] = coord

    # 校验 mapping
    for old_key, new_key in mapping_old_to_new.items():
        if old_key not in old_dict:
            logger.print(f"[ERROR] Old residue not found: {old_key}")
            return False

        if new_key not in new_dict:
            logger.print(f"[ERROR] New residue not found: {new_key}")
            return False

        old_coord = old_dict[old_key]
        new_coord = new_dict[new_key]

        dist = math.sqrt(
            (old_coord[0] - new_coord[0]) ** 2 +
            (old_coord[1] - new_coord[1]) ** 2 +
            (old_coord[2] - new_coord[2]) ** 2
        )

        if dist > tol:
            logger.print(
                f"[ERROR] CA coordinate changed after cleaning: "
                f"{old_key} -> {new_key}, distance = {dist}"
            )
            return False

    return True



def generate_clean_report(structure: Structure, cleaned_structure: Structure, mapping_old_to_new: Dict[Tuple[int, str, str], Tuple[int, str, str]], stats: Dict[str, int], logger: Logger) -> dict | None:

    amino_acid_mapping_old_2_new: list[dict[str, dict]] = []

    old_chain = get_single_chain(structure, logger)
    if old_chain is None:
        return None

    new_chain = get_single_chain(cleaned_structure, logger)
    if new_chain is None:
        return None

    old_residue_dict = {}
    for res in old_chain.get_residues():
        hetflag, resseq, icode = res.id
        if hetflag != " ":
            continue
        resname = res.get_resname().strip()
        old_residue_dict[(int(resseq), resname, str(icode))] = res

    new_residue_dict = {}
    for res in new_chain.get_residues():
        hetflag, resseq, icode = res.id
        if hetflag != " ":
            continue
        resname = res.get_resname().strip()
        new_residue_dict[(int(resseq), resname, str(icode))] = res

    for old_key, new_value in mapping_old_to_new.items():
        old_resseq, old_resname, old_icode = old_key
        new_resseq, new_resname, new_icode = new_value

        old_res = old_residue_dict.get((old_resseq, old_resname, old_icode))
        if old_res is None:
            logger.print(f"[ERROR] Old residue not found: {old_key}")
            return None

        new_res = new_residue_dict.get((new_resseq, new_resname, new_icode))
        if new_res is None:
            logger.print(f"[ERROR] New residue not found: {new_value}")
            return None

        old_h_count = sum(
            1 for atom in old_res.get_atoms()
            if atom.element and atom.element.strip().upper() == "H"
        )

        new_h_count = sum(
            1 for atom in new_res.get_atoms()
            if atom.element and atom.element.strip().upper() == "H"
        )

        amino_acid_mapping_old_2_new.append({
            "old_residue": {
                "aa_id": old_resseq,
                "aa_name": normalize_aa_name_to_one_letter(old_resname),
                "hydrogen_atom_count": old_h_count,
            },
            "new_residue": {
                "aa_id": new_resseq,
                "aa_name": normalize_aa_name_to_one_letter(new_resname),
                "hydrogen_atom_count": new_h_count,
            }
        })

    return {
        "output_type": "enzywizard_clean",
        "amino_acid_mapping_old_to_new": amino_acid_mapping_old_2_new,
        "clean_statistics": stats,
    }

def check_cleaned_structure(struct: Structure, logger: Logger) -> bool:
    model_count = len(list(struct.get_models()))
    if model_count != 1:
        logger.print("[ERROR] Structure must contain exactly one model. Please run 'enzywizard clean' first.")
        return False

    model = next(struct.get_models())
    chains = list(model.get_chains())

    if len(chains) != 1:
        logger.print("[ERROR] Structure must contain exactly one chain. Please run 'enzywizard clean' first.")
        return False

    chain = chains[0]
    if chain.id != "A":
        logger.print("[ERROR] Cleaned structure must use chain ID 'A'. Please run 'enzywizard clean' first.")
        return False

    expected_resseq = 1

    for res in chain.get_residues():
        hetflag, resseq, icode = res.id

        if str(hetflag).strip():
            logger.print(f"[ERROR] Non-protein or hetero residue detected at residue {resseq}{icode}. Please run 'enzywizard clean' first.")
            return False

        if str(icode).strip():
            logger.print(f"[ERROR] Insertion code detected at residue {resseq}{icode}. Please run 'enzywizard clean' first.")
            return False

        resname = res.get_resname().strip()
        resname_std = standardize_resname(resname)
        if resname != resname_std:
            logger.print(f"[ERROR] Residue name '{resname}' at residue {resseq} is not standardized. Please run 'enzywizard clean' first.")
            return False

        if resname not in AA3_STANDARD:
            logger.print(f"[ERROR] Non-standard residue '{resname}' detected at residue {resseq}. Please run 'enzywizard clean' first.")
            return False

        if int(resseq) != expected_resseq:
            logger.print(
                f"[ERROR] Residue numbering is not continuous: expected {expected_resseq}, got {resseq}. Please run 'enzywizard clean' first."
            )
            return False

        atoms_by_name: Dict[str, List[Atom]] = {}
        for atom in res.get_atoms():
            if is_hydrogen_atom(atom):
                continue

            atom_name = normalize_atom_name(atom.get_name())
            atoms_by_name.setdefault(atom_name, []).append(atom)

        for atom_name in BACKBONE_REQUIRED_ATOMS:
            if atom_name not in atoms_by_name:
                logger.print(
                    f"[ERROR] Missing backbone atom '{atom_name}' at residue {resseq}. Please run 'enzywizard clean' first.")
                return False

        for atom_name in BACKBONE_REQUIRED_ATOMS:
            chosen = choose_atom_altloc(atoms_by_name[atom_name])
            occ = chosen.get_occupancy()
            if occ is not None and occ < 0:
                logger.print(
                    f"[ERROR] Backbone atom '{atom_name}' at residue {resseq} has invalid occupancy {occ}. Please run 'enzywizard clean' first."
                )
                return False

        expected_heavy_atom_set = AA3_EXPECTED_HEAVY_ATOM_SET.get(resname)
        if expected_heavy_atom_set is None:
            logger.print(f"[ERROR] No expected heavy atom definition for residue '{resname}' at residue {resseq}.")
            return False

        actual_heavy_atom_set = set(atoms_by_name.keys())

        if not expected_heavy_atom_set.issubset(actual_heavy_atom_set):
            missing_atom_name_list = sorted(expected_heavy_atom_set - actual_heavy_atom_set)
            logger.print(
                f"[ERROR] Missing heavy atoms at residue {resseq} ({resname}): {missing_atom_name_list}. "
                f"Please run 'enzywizard clean' first."
            )
            return False

        allowed_heavy_atom_set = AA3_ALLOWED_HEAVY_ATOM_SET_WITH_OXT[resname]
        if not actual_heavy_atom_set.issubset(allowed_heavy_atom_set):
            unexpected_atom_name_list = sorted(actual_heavy_atom_set - allowed_heavy_atom_set)
            logger.print(
                f"[ERROR] Unexpected heavy atoms at residue {resseq} ({resname}): {unexpected_atom_name_list}. "
                f"Please run 'enzywizard clean' first."
            )
            return False

        for atom_name in actual_heavy_atom_set:
            chosen = choose_atom_altloc(atoms_by_name[atom_name])
            occ = chosen.get_occupancy()
            if occ is not None and occ < 0:
                logger.print(
                    f"[ERROR] Heavy atom '{atom_name}' at residue {resseq} has invalid occupancy {occ}. Please run 'enzywizard clean' first."
                )
                return False

        expected_resseq += 1

    return True