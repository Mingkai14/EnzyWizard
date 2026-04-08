from __future__ import annotations
from typing import List
from ..utils.logging_utils import Logger

def moving_average(values: List[float], window_size: int, logger: Logger) -> List[float] | None:
    if window_size <= 0:
        logger.print("[ERROR] window_size must be positive")
        return None

    n = len(values)
    if n == 0:
        logger.print("[ERROR] Empty input for moving_average")
        return None

    half = window_size // 2
    result: List[float] = []

    for i in range(n):
        left = max(0, i - half)
        right = min(n, i + half + 1)
        sub = values[left:right]
        result.append(sum(sub) / len(sub))

    return result
