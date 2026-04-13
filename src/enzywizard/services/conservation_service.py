from __future__ import annotations
from ..utils.logging_utils import Logger
from pathlib import Path
from ..utils.IO_utils import file_exists, get_stem, check_filename_length, load_msa, load_fasta, write_msa, write_hmm, write_json_from_dict_inline_leaf_lists
from ..utils.sequence_utils import check_msa, clean_msa_to_sto
from ..algorithms.conservation_algorithms import compute_conservation_scores, generate_conservation_report
from ..utils.common_utils import get_optimized_filename

def run_conservation_service(input_fasta: str | Path, input_msa: str | Path, output_dir: str | Path) ->bool:
    # ---- logger ----
    logger = Logger(output_dir)
    logger.print(f"[INFO] Conservation processing started: {input_fasta} {input_msa}")

    # ---- check input ----
    input_fasta = Path(input_fasta)
    input_msa = Path(input_msa)
    output_dir = Path(output_dir)

    if not file_exists(input_fasta):
        logger.print(f"[ERROR] Input not found: {input_fasta}")
        return False

    if not file_exists(input_msa):
        logger.print(f"[ERROR] Input not found: {input_msa}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- get name ----
    protein_name = get_stem(input_fasta)
    if not check_filename_length(protein_name,logger):
        return False
    logger.print(f"[INFO] Protein name resolved: {protein_name}")

    msa_name = get_stem(input_msa)
    if not check_filename_length(msa_name,logger):
        return False
    logger.print(f"[INFO] MSA file name resolved: {msa_name}")

    # ---- load files ----
    sequence_dict = load_fasta(input_fasta,logger)
    if sequence_dict is None:
        return False

    logger.print("[INFO] FASTA file loaded")

    msa_list = load_msa(input_msa,logger)
    if msa_list is None:
        return False

    logger.print("[INFO] MSA file loaded")

    #---- check msa ----
    if not check_msa(input_msa,sequence_dict, msa_list, logger):
        return False
    logger.print(f"[INFO] MSA checked")

    #---- clean msa ----
    cleaned_msa_list=clean_msa_to_sto(msa_list,logger)
    if cleaned_msa_list is None:
        return False
    logger.print(f"[INFO] MSA identically cleaned to STO format")

    #---- write cleaned msa ----
    cleaned_msa_path=output_dir / get_optimized_filename(f'cleaned_{msa_name}.sto')
    if not write_msa(cleaned_msa_list,cleaned_msa_path,logger):
        return False
    logger.print(f"[INFO] Cleaned MSA STO file saved: {str(cleaned_msa_path)}")

    #---- write hmm file ----
    hmm_path = output_dir / get_optimized_filename(f'hmm_profile_{msa_name}.hmm')
    if not write_hmm(cleaned_msa_path,hmm_path,logger):
        return False
    logger.print(f"[INFO] HMM Profile file saved: {str(hmm_path)}")

    #---- run algorithm ----
    logger.print("[INFO] Calculating conservation scores started")
    conservation_scores=compute_conservation_scores(hmm_path,sequence_dict,logger)
    if conservation_scores is None:
        return False

    json_report_path = output_dir / get_optimized_filename(f"conservation_report_{protein_name}.json")

    report=generate_conservation_report(conservation_scores)
    write_json_from_dict_inline_leaf_lists(report,json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Conservation processing finished")

    return True









