from __future__ import annotations
from Bio.PDB.Atom import Atom
from typing import Dict, List, Tuple, Optional, Union
from ..resources.aa_resources import AA3_STANDARD, modres

def standardize_resname(resname: str) -> str:
    resname = resname.strip()
    if resname in AA3_STANDARD:
        return resname
    if resname in modres:
        return modres[resname]
    return resname

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
    return Atom(
        name=atom.get_name(),
        coord=coord,
        bfactor=atom.get_bfactor(),
        occupancy=atom.get_occupancy(),
        altloc=" ",  # remove altloc
        fullname=atom.get_fullname(),
        serial_number=atom.get_serial_number(),
        element=atom.element,
    )