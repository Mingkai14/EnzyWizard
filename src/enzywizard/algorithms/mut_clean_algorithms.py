from typing import Dict, Tuple, Any
from ..utils.logging_utils import Logger
from ..utils.mut_clean_utils import get_muts_from_aas
from Bio.PDB.Structure import Structure
from ..utils.structure_utils import get_single_chain
from ..utils.sequence_utils import normalize_aa_name_to_one_letter

def get_cleaned_amino_acid_substitution(
    wt_mapping_old_to_new: Dict[Tuple[int, str, str], Tuple[int, str, str]],
    mut_mapping_old_to_new: Dict[Tuple[int, str, str], Tuple[int, str, str]],
    mutation: str,
    logger: Logger,
) -> str | None:


    muts_list = get_muts_from_aas(mutation)

    # 1. mapping 长度一致
    if len(wt_mapping_old_to_new) != len(mut_mapping_old_to_new):
        logger.print(
            f"[ERROR] Wild-type and mutant mappings have different lengths: {len(wt_mapping_old_to_new)} vs {len(mut_mapping_old_to_new)}."
        )
        return None

    # 建索引
    wt_pos_to_new = {}
    mut_pos_to_new = {}

    for (old_resseq, resname_old, old_icode), (new_resseq, resname_std, new_icode) in wt_mapping_old_to_new.items():
        wt_pos_to_new[old_resseq] = new_resseq

    for (old_resseq, resname_old, old_icode), (new_resseq, resname_std, new_icode) in mut_mapping_old_to_new.items():
        mut_pos_to_new[old_resseq] = new_resseq

    cleaned_muts = []

    for wt_aa,old_pos,mut_aa in muts_list:
        # 2. 每个单点突的位置都要在两个 mapping 中找到
        if old_pos not in wt_pos_to_new:
            logger.print(f"[ERROR] Mutation position {old_pos} not found in wild-type mapping.")
            return None

        if old_pos not in mut_pos_to_new:
            logger.print(f"[ERROR] Mutation position {old_pos} not found in mutant mapping.")
            return None

        wt_new_resseq = wt_pos_to_new[old_pos]
        mut_new_resseq = mut_pos_to_new[old_pos]

        # 3. WT 和 MUT 对应的新序号必须一致
        if wt_new_resseq != mut_new_resseq:
            logger.print(
                f"[ERROR] Inconsistent cleaned residue position for original mutation position {old_pos}: wild-type -> {wt_new_resseq}, mutant -> {mut_new_resseq}."
            )
            return None

        cleaned_muts.append(f"{wt_aa}{wt_new_resseq}{mut_aa}")

    cleaned_mutation = ",".join(cleaned_muts)
    return cleaned_mutation


def generate_mutclean_report(
    old_aas: str,
    cleaned_aas: str,
    wt_structure: Structure,
    wt_cleaned_structure: Structure,
    mut_structure: Structure,
    mut_cleaned_structure: Structure,
    wt_mapping_old_to_new: Dict[Tuple[int, str, str], Tuple[int, str, str]],
    wt_stats: Dict[str, int],
    mut_mapping_old_to_new: Dict[Tuple[int, str, str], Tuple[int, str, str]],
    mut_stats: Dict[str, int],
    logger: Logger
) -> dict | None:
    wt_amino_acid_mapping_old_2_new: list[dict[str, dict]] = []
    mut_amino_acid_mapping_old_2_new: list[dict[str, dict]] = []

    wt_old_chain = get_single_chain(wt_structure, logger)
    if wt_old_chain is None:
        return None

    wt_new_chain = get_single_chain(wt_cleaned_structure, logger)
    if wt_new_chain is None:
        return None

    mut_old_chain = get_single_chain(mut_structure, logger)
    if mut_old_chain is None:
        return None

    mut_new_chain = get_single_chain(mut_cleaned_structure, logger)
    if mut_new_chain is None:
        return None

    wt_old_residue_dict = {}
    for res in wt_old_chain.get_residues():
        hetflag, resseq, icode = res.id
        if hetflag != " ":
            continue
        resname = res.get_resname().strip()
        wt_old_residue_dict[(int(resseq), resname, str(icode))] = res

    wt_new_residue_dict = {}
    for res in wt_new_chain.get_residues():
        hetflag, resseq, icode = res.id
        if hetflag != " ":
            continue
        resname = res.get_resname().strip()
        wt_new_residue_dict[(int(resseq), resname, str(icode))] = res

    mut_old_residue_dict = {}
    for res in mut_old_chain.get_residues():
        hetflag, resseq, icode = res.id
        if hetflag != " ":
            continue
        resname = res.get_resname().strip()
        mut_old_residue_dict[(int(resseq), resname, str(icode))] = res

    mut_new_residue_dict = {}
    for res in mut_new_chain.get_residues():
        hetflag, resseq, icode = res.id
        if hetflag != " ":
            continue
        resname = res.get_resname().strip()
        mut_new_residue_dict[(int(resseq), resname, str(icode))] = res

    for old_key, new_value in wt_mapping_old_to_new.items():
        old_resseq, old_resname, old_icode = old_key
        new_resseq, new_resname, new_icode = new_value

        old_res = wt_old_residue_dict.get((old_resseq, old_resname, old_icode))
        if old_res is None:
            logger.print(f"[ERROR] WT old residue not found: {old_key}")
            return None

        new_res = wt_new_residue_dict.get((new_resseq, new_resname, new_icode))
        if new_res is None:
            logger.print(f"[ERROR] WT new residue not found: {new_value}")
            return None

        old_h_count = sum(
            1 for atom in old_res.get_atoms()
            if atom.element and atom.element.strip().upper() == "H"
        )

        new_h_count = sum(
            1 for atom in new_res.get_atoms()
            if atom.element and atom.element.strip().upper() == "H"
        )

        wt_amino_acid_mapping_old_2_new.append({
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

    for old_key, new_value in mut_mapping_old_to_new.items():
        old_resseq, old_resname, old_icode = old_key
        new_resseq, new_resname, new_icode = new_value

        old_res = mut_old_residue_dict.get((old_resseq, old_resname, old_icode))
        if old_res is None:
            logger.print(f"[ERROR] MUT old residue not found: {old_key}")
            return None

        new_res = mut_new_residue_dict.get((new_resseq, new_resname, new_icode))
        if new_res is None:
            logger.print(f"[ERROR] MUT new residue not found: {new_value}")
            return None

        old_h_count = sum(
            1 for atom in old_res.get_atoms()
            if atom.element and atom.element.strip().upper() == "H"
        )

        new_h_count = sum(
            1 for atom in new_res.get_atoms()
            if atom.element and atom.element.strip().upper() == "H"
        )

        mut_amino_acid_mapping_old_2_new.append({
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
        "output_type": "enzywizard_mut_clean",
        "amino_acid_substitution": old_aas,
        "cleaned_amino_acid_substitution": cleaned_aas,
        "wt_amino_acid_mapping_old_to_new": wt_amino_acid_mapping_old_2_new,
        "wt_clean_statistics": wt_stats,
        "mut_amino_acid_mapping_old_to_new": mut_amino_acid_mapping_old_2_new,
        "mut_clean_statistics": mut_stats,
    }