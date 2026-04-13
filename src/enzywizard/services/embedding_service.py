from __future__ import annotations
from pathlib import Path
from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists,get_stem, check_filename_length, load_fasta, write_json_from_dict_inline_leaf_lists
from ..algorithms.embedding_algorithms import generate_embedding, generate_embedding_report
from ..utils.common_utils import get_optimized_filename

def run_embedding_service(input_fasta: str | Path, output_dir: str | Path, model_name: str = "esm2_t6_8M_UR50D") -> bool:
    # ---- logger ----
    logger = Logger(output_dir)
    logger.print(f"[INFO] Embedding processing started: {input_fasta}")

    # ---- check input ----
    input_path = Path(input_fasta)
    output_dir = Path(output_dir)

    if not file_exists(input_path):
        logger.print(f"[ERROR] Input not found: {input_path}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- get name ----
    name = get_stem(input_path)
    if not check_filename_length(name, logger):
        return False
    logger.print(f"[INFO] Protein name resolved: {name}")

    # ---- load sequence ----
    sequence_dict = load_fasta(input_fasta,logger)
    if sequence_dict is None:
        return False

    logger.print("[INFO] FASTA file loaded")

    # ---- run algorithm ----
    logger.print("[INFO] Generating protein embeddings started")
    embeddings=generate_embedding(sequence_dict,logger,model_name=model_name)
    if embeddings is None:
        return False

    # ---- generate report ----
    report = generate_embedding_report(embeddings)

    # ---- write output ----
    json_report_path = output_dir / get_optimized_filename(f"embedding_report_{name}.json")
    write_json_from_dict_inline_leaf_lists(report, json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Embedding processing finished")

    return True

