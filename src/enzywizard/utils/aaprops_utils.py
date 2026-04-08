from Bio.PDB.DSSP import DSSP
from typing import Tuple
from ..resources.aa_resources import DSSP_8STATE
from ..utils.logging_utils import Logger

def get_dssp_fields(dssp: DSSP, dssp_key: Tuple[str, Tuple[str, int, str]], logger: Logger) -> Tuple[str, float, float, float] | None:

    if dssp_key not in dssp:
        logger.print(f"[ERROR] DSSP key {dssp_key} not found.")
        return None

    dssp_res = dssp[dssp_key]

    ss = dssp_res[2]
    rsa = dssp_res[3]
    phi = dssp_res[4]
    psi = dssp_res[5]

    ss = (ss or "-").strip()
    if ss == "" or ss not in DSSP_8STATE:
        ss = "-"

    try:
        return ss, float(rsa), float(phi), float(psi)
    except (TypeError, ValueError):
        logger.print(f"[ERROR] Invalid DSSP numeric values at {dssp_key}: rsa={rsa}, phi={phi}, psi={psi}")
        return None