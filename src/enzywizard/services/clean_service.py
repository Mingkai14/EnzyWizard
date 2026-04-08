from __future__ import annotations
from ..utils.logging_utils import Logger
from pathlib import Path
from ..utils.IO_utils import file_exists, get_stem, check_filename_length, load_protein_structure, write_cif ,write_pdb, write_json_from_dict_inline_leaf_lists,structure_to_pdbfile,modeller_to_structure, write_fasta
from ..algorithms.clean_algorithms import clean_structure_to_single_chain_A, generate_clean_report, check_cleaned_structure, add_hydrogens_to_pdbfile, validate_clean_mapping_coordinates

def run_clean_service(input_path: str | Path, output_dir: str | Path, add_H: bool = True, pH: float = 7.0, force_field_file: str = "charmm36.xml") ->bool:
    # ---- logger ----
    logger = Logger(output_dir)
    logger.print(f"[INFO] Clean processing started: {input_path}")

    # ---- check input ----
    input_path = Path(input_path)
    output_dir = Path(output_dir)

    if not file_exists(input_path):
        logger.print(f"[ERROR] Input not found: {input_path}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- get name ----
    name = get_stem(input_path)
    if not check_filename_length(name,logger):
        return False
    logger.print(f"[INFO] Protein name resolved: {name}")

    # ---- load structure ----
    structure = load_protein_structure(input_path,name,logger)
    if structure is None:
        logger.print(f"[ERROR] Failed to load structure: {input_path}")
        return False

    logger.print("[INFO] Structure loaded")

    # ---- run algorithm ----
    logger.print("[INFO] Cleaning started")
    clean_result = clean_structure_to_single_chain_A(structure, logger)
    if clean_result is None:
        return False
    cleaned_structure, mapping_old_to_new, stats= clean_result

    if add_H:
        temp_pdbfile=structure_to_pdbfile(cleaned_structure,logger,protein_name=name)
        if temp_pdbfile is None:
            return False
        temp_modeller=add_hydrogens_to_pdbfile(temp_pdbfile,logger,pH=pH,force_field_file=force_field_file)
        if temp_modeller is None:
            return False
        cleaned_structure=modeller_to_structure(temp_modeller,logger,protein_name=name)
        if cleaned_structure is None:
            return False
        logger.print(f"[INFO] Hydrogens added")

    #---- check structure ----
    if not check_cleaned_structure(cleaned_structure, logger):
        return False
    if not validate_clean_mapping_coordinates(structure,cleaned_structure,mapping_old_to_new,logger):
        return False

    logger.print(f"[INFO] Structure checked")

    # ---- save results ----

    cleaned_cif_path = output_dir / f"cleaned_{name}.cif"
    cleaned_pdb_path = output_dir / f"cleaned_{name}.pdb"
    cleaned_fasta_path = output_dir / f"cleaned_{name}.fasta"
    json_report_path = output_dir / f"clean_report_{name}.json"

    write_cif(cleaned_structure,cleaned_cif_path)
    logger.print(f"[INFO] Cleaned CIF saved: {cleaned_cif_path}")

    write_pdb(cleaned_structure,cleaned_pdb_path)
    logger.print(f"[INFO] Cleaned PDB saved: {cleaned_pdb_path}")

    if write_fasta(cleaned_structure,cleaned_fasta_path, logger):
        logger.print(f"[INFO] Cleaned FASTA saved: {cleaned_fasta_path}")
    else:
        return False

    report=generate_clean_report(structure,cleaned_structure,mapping_old_to_new,stats,logger)
    if report is None:
        return False
    write_json_from_dict_inline_leaf_lists(report,json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Clean processing finished")

    return True
