import re
from ..utils.logging_utils import Logger
from typing import List, Tuple

def check_amino_acid_substitution(mutation: str,wt_length: int,mut_length: int,logger: Logger) -> bool:

    if not isinstance(mutation, str) or not mutation.strip():
        logger.print("[ERROR] Empty amino acid substitution.")
        return False

    muts = mutation.split(",")
    seen_positions = set()
    for single_mut in muts:
        single_mut = single_mut.strip()

        m = re.fullmatch(
            r"([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])",
            single_mut,
            flags=re.I,
        )
        if not m:
            logger.print(f"[ERROR] Invalid amino acid substitution format: {single_mut}.")
            return False

        wt_aa = m.group(1).upper()
        pos_str = m.group(2)
        mut_aa = m.group(3).upper()

        if wt_aa == mut_aa:
            logger.print(f"[ERROR] No amino acid mutated: {single_mut}.")
            return False

        pos = int(pos_str)

        if pos in seen_positions:
            logger.print(f"[ERROR] Duplicate mutation position detected: {pos}.")
            return False
        seen_positions.add(pos)

        if pos < 1 or pos > wt_length:
            logger.print(f"[ERROR] Mutation position out of wild-type range: {single_mut} (position={pos}, wt_length={wt_length})")
            return False

        if pos < 1 or pos > mut_length:
            logger.print(f"[ERROR] Mutation position out of mutant range: {single_mut} (position={pos}, mut_length={mut_length})")
            return False
    return True

def get_muts_from_aas(mutation: str)->List[Tuple[str, int, str]]:
    muts_list: List[Tuple[str, int, str]]=[]
    muts = mutation.split(",")

    for single_mut in muts:
        single_mut = single_mut.strip()

        m = re.fullmatch(
            r"([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])",
            single_mut,
            flags=re.I,
        )

        wt_aa = m.group(1).upper()
        pos_str = m.group(2)
        mut_aa = m.group(3).upper()
        muts_list.append((wt_aa,int(pos_str),mut_aa))
    return muts_list
