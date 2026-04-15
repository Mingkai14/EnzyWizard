from __future__ import annotations
from Bio.PDB.Atom import Atom
from typing import Dict, List, Tuple, Optional, Union
from Bio.PDB.Residue import Residue
from ..resources.aa_resources import AA3_STANDARD, modres

def standardize_resname(resname: str) -> str:
    resname = resname.strip().upper()
    if resname in AA3_STANDARD:
        return resname
    if resname in modres:
        return modres[resname].strip().upper()
    return resname

def normalize_atom_name(atom_name: str) -> str:
    return atom_name.strip().upper()

def is_hydrogen_atom(atom: Atom) -> bool:
    element = atom.element.strip().upper() if atom.element is not None else ""
    name = normalize_atom_name(atom.get_name())

    if element == "H":
        return True
    if name.startswith("H"):
        return True
    return False

def get_residue_heavy_atom_name_set(res: Residue) -> set[str]:
    heavy_atom_name_set: set[str] = set()

    for atom in res.get_atoms():
        if is_hydrogen_atom(atom):
            continue
        heavy_atom_name_set.add(normalize_atom_name(atom.get_name()))

    return heavy_atom_name_set

def choose_atom_altloc(atom_list: List[Atom]) -> Atom:
    # Prefer blank altloc
    for a in atom_list:
        if (a.get_altloc() or " ").strip() == "":
            return a

    # Else pick highest occupancy
    best = atom_list[0]
    best_occ = best.get_occupancy()
    best_occ = best_occ if best_occ is not None else -1.0
    for a in atom_list[1:]:
        occ = a.get_occupancy()
        occ = occ if occ is not None else -1.0
        if occ > best_occ:
            best, best_occ = a, occ
    return best

def clone_atom(atom: Atom, *, new_coord=None) -> Atom:
    coord = new_coord if new_coord is not None else atom.get_coord()
    normalized_name = normalize_atom_name(atom.get_name())

    return Atom(
        name=normalized_name,
        coord=coord,
        bfactor=atom.get_bfactor(),
        occupancy=atom.get_occupancy(),
        altloc=" ",
        fullname=f"{normalized_name:>4}",
        serial_number=atom.get_serial_number(),
        element=atom.element.strip().upper() if atom.element is not None else atom.element,
    )