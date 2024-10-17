from typing import Any, Optional
import os
import pickle as pkl

from autogrow.utils.logging import log_info

def _get_cache_filename(id: str, cache_dir: Optional[str] = None):
    return (
        os.path.join(cache_dir, f"{id}_results_cache.pkl")
        if cache_dir is not None
        else ""
    )

def load_from_cache(id: str, cache_dir: str) -> Optional[Any]:
    cache_filename = _get_cache_filename(id, cache_dir=cache_dir)

    if not os.path.exists(cache_filename):
        return None

    log_info(f"Loading previous {id} results from cache: {cache_filename}")
    return pkl.load(open(cache_filename, "rb"))

def save_to_cache(id: str, cache_dir: str, data: Any) -> None:
    cache_filename = _get_cache_filename(id, cache_dir=cache_dir)
    log_info(f"Saving {id} results to cache: {cache_filename}")
    pkl.dump(data, open(cache_filename, "wb"))