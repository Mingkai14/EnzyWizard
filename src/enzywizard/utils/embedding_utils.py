from __future__ import annotations
import torch
import esm
from typing import Tuple
from functools import lru_cache

@lru_cache(maxsize=4)
def load_esm2(model_name: str = "esm2_t6_8M_UR50D") -> Tuple[torch.nn.Module, esm.data.Alphabet]:
    """
    Load ESM2 model + alphabet once (cached).
    """
    loader = getattr(esm.pretrained, model_name, None)
    if loader is None:
        raise ValueError(f"Unknown ESM2 model name: {model_name}")

    model, alphabet = loader()
    model.eval()
    return model, alphabet