from __future__ import annotations
from pathlib import Path
from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists, get_stem, check_filename_length, load_protein_structure, write_cif, write_pdb, write_json_from_dict_inline_leaf_lists, structure_to_pdbfile, modeller_to_structure, write_fasta
from ..utils.mut_clean_utils import check_amino_acid_substitution
from ..utils.structure_utils import get_single_chain,get_chain_length
from ..algorithms.clean_algorithms import clean_structure_to_single_chain_A, add_hydrogens_to_pdbfile, check_cleaned_structure, validate_clean_mapping_coordinates
from ..algorithms.mut_clean_algorithms import get_cleaned_amino_acid_substitution, generate_mutclean_report
from ..utils.common_utils import get_optimized_filename


def run_mutclean_service(wt_input_path: str | Path, mut_input_path: str | Path, mutation:str, output_dir: str | Path, add_H: bool = True, pH: float = 7.0, force_field_file: str = "charmm36.xml")->bool:
    # ---- logger ----
    logger = Logger(output_dir)
    logger.print(f"[INFO] Mutclean processing started: {wt_input_path} {mut_input_path} {mutation}")

    # ---- check input ----
    wt_input_path=Path(wt_input_path)
    mut_input_path=Path(mut_input_path)
    output_dir = Path(output_dir)

    if not file_exists(wt_input_path):
        logger.print(f"[ERROR] Input not found: {wt_input_path}")
        return False

    if not file_exists(mut_input_path):
        logger.print(f"[ERROR] Input not found: {mut_input_path}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- get name ----
    wt_name = get_stem(wt_input_path)
    if not check_filename_length(wt_name,logger):
        return False
    mut_name = get_stem(mut_input_path)
    if not check_filename_length(mut_name,logger):
        return False
    if wt_name==mut_name:
        logger.print(f"[ERROR] Wild-type and mutant names are the same: {wt_name} {mut_name}")
        return False
    logger.print(f"[INFO] Protein names resolved: {wt_name} {mut_name}")

    # ---- load structure ----
    wt_structure = load_protein_structure(wt_input_path,wt_name,logger)
    mut_structure = load_protein_structure(mut_input_path,mut_name,logger)
    if wt_structure is None:
        logger.print(f"[ERROR] Failed to load wild-type structure: {wt_input_path}")
        return False
    if mut_structure is None:
        logger.print(f"[ERROR] Failed to load mutant structure: {mut_input_path}")
        return False

    logger.print("[INFO] Structures loaded")

    # ---- check amino acid substitution ----
    wt_chain = get_single_chain(wt_structure, logger)
    mut_chain = get_single_chain(mut_structure, logger)
    if wt_chain is None or mut_chain is None:
        return False

    wt_seq_length=get_chain_length(wt_chain,logger)
    mut_seq_length=get_chain_length(mut_chain,logger)
    if wt_seq_length is None or mut_seq_length is None:
        return False
    if not check_amino_acid_substitution(mutation,wt_seq_length,mut_seq_length,logger):
        return False

    # ---- check length ----
    if wt_seq_length!=mut_seq_length:
        logger.print(f"[ERROR] Wild-type and mutant sequence lengths are not equal: {wt_seq_length} {mut_seq_length}")
        return False

    # ---- run algorithm ----
    logger.print("[INFO] Cleaning wild-type and mutant structures started")
    clean_result = clean_structure_to_single_chain_A(wt_structure, logger)
    if clean_result is None:
        return False
    wt_cleaned_structure, wt_mapping_old_to_new, wt_stats= clean_result


    clean_result = clean_structure_to_single_chain_A(mut_structure, logger)
    if clean_result is None:
        return False
    mut_cleaned_structure, mut_mapping_old_to_new, mut_stats= clean_result

    cleaned_mutation=get_cleaned_amino_acid_substitution(wt_mapping_old_to_new,mut_mapping_old_to_new,mutation,logger)
    if cleaned_mutation is None:
        return False

    if add_H:
        temp_pdbfile=structure_to_pdbfile(wt_cleaned_structure,logger,protein_name=wt_name)
        if temp_pdbfile is None:
            return False
        temp_modeller=add_hydrogens_to_pdbfile(temp_pdbfile,logger,pH=pH,force_field_file=force_field_file)
        if temp_modeller is None:
            return False
        wt_cleaned_structure=modeller_to_structure(temp_modeller,logger,protein_name=wt_name)
        if wt_cleaned_structure is None:
            return False

        temp_pdbfile=structure_to_pdbfile(mut_cleaned_structure,logger,protein_name=mut_name)
        if temp_pdbfile is None:
            return False
        temp_modeller=add_hydrogens_to_pdbfile(temp_pdbfile,logger,pH=pH,force_field_file=force_field_file)
        if temp_modeller is None:
            return False
        mut_cleaned_structure=modeller_to_structure(temp_modeller,logger,protein_name=mut_name)
        if mut_cleaned_structure is None:
            return False

        logger.print(f"[INFO] Hydrogens added")

    # ---- post check ----
    wt_cleaned_chain = get_single_chain(wt_cleaned_structure, logger)
    mut_cleaned_chain = get_single_chain(mut_cleaned_structure, logger)
    if wt_cleaned_chain is None or mut_cleaned_chain is None:
        return False

    wt_cleaned_length=get_chain_length(wt_cleaned_chain,logger)
    mut_cleaned_length=get_chain_length(mut_cleaned_chain,logger)
    if wt_cleaned_length!=mut_cleaned_length:
        logger.print(f"[ERROR] Cleaned wild-type and mutant sequence lengths are not equal: {wt_cleaned_length} {mut_cleaned_length}")
        return False

    if not check_cleaned_structure(wt_cleaned_structure, logger):
        return False
    if not check_cleaned_structure(mut_cleaned_structure, logger):
        return False
    if not validate_clean_mapping_coordinates(wt_structure,wt_cleaned_structure,wt_mapping_old_to_new,logger):
        return False
    if not validate_clean_mapping_coordinates(mut_structure,mut_cleaned_structure,mut_mapping_old_to_new,logger):
        return False

    # ---- save results ----

    cleaned_wt_cif_path = output_dir / get_optimized_filename(f"cleaned_{wt_name}.cif")
    cleaned_wt_pdb_path = output_dir / get_optimized_filename(f"cleaned_{wt_name}.pdb")
    cleaned_wt_fasta_path = output_dir / get_optimized_filename(f"cleaned_{wt_name}.fasta")

    cleaned_mut_cif_path = output_dir / get_optimized_filename(f"cleaned_{mut_name}.cif")
    cleaned_mut_pdb_path = output_dir / get_optimized_filename(f"cleaned_{mut_name}.pdb")
    cleaned_mut_fasta_path = output_dir / get_optimized_filename(f"cleaned_{mut_name}.fasta")

    json_report_path = output_dir / get_optimized_filename(f"mut_clean_report_{wt_name}_to_{mut_name}.json")

    write_cif(wt_cleaned_structure,cleaned_wt_cif_path)
    logger.print(f"[INFO] Cleaned wild-type CIF saved: {cleaned_wt_cif_path}")

    write_pdb(wt_cleaned_structure,cleaned_wt_pdb_path)
    logger.print(f"[INFO] Cleaned wild-type PDB saved: {cleaned_wt_pdb_path}")

    if write_fasta(wt_cleaned_structure,cleaned_wt_fasta_path, logger):
        logger.print(f"[INFO] Cleaned wild-type FASTA saved: {cleaned_wt_fasta_path}")
    else:
        return False

    write_cif(mut_cleaned_structure,cleaned_mut_cif_path)
    logger.print(f"[INFO] Cleaned mutant CIF saved: {cleaned_mut_cif_path}")

    write_pdb(mut_cleaned_structure,cleaned_mut_pdb_path)
    logger.print(f"[INFO] Cleaned mutant PDB saved: {cleaned_mut_pdb_path}")

    if write_fasta(mut_cleaned_structure,cleaned_mut_fasta_path, logger):
        logger.print(f"[INFO] Cleaned mutant FASTA saved: {cleaned_mut_fasta_path}")
    else:
        return False

    report=generate_mutclean_report(mutation,cleaned_mutation,wt_structure,wt_cleaned_structure,mut_structure,mut_cleaned_structure,wt_mapping_old_to_new,wt_stats,mut_mapping_old_to_new,mut_stats,logger)
    if report is None:
        return False
    write_json_from_dict_inline_leaf_lists(report,json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Mutclean processing finished")

    return True
