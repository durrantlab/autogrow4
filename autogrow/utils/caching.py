from typing import Any, Optional
import os
import pickle as pkl
from autogrow.utils.logging import log_debug, log_info

class CacheManager:
    def __init__(self, label: str, cache_dir: str):
        self.label = label
        self.cache_dir = cache_dir
        self.exists = False
        self.data: Any = None

    def _get_cache_filename(self) -> str:
        return os.path.join(self.cache_dir, f"{self.label}_results_cache.pkl")

    def _load_from_cache(self) -> Optional[Any]:
        cache_filename = self._get_cache_filename()
        self.exists = os.path.exists(cache_filename)
        if not self.exists:
            return None
        with open(cache_filename, "rb") as f:
            c = pkl.load(f)
            log_info(f"Loaded previous {self.label} results from cache: {cache_filename}")
            return c

    def _save_to_cache(self) -> None:
        if not self.exists and self.data is not None:
            cache_filename = self._get_cache_filename()
            # TODO: What if dir of cache_filename doesn't exist? Happened once. 
            with open(cache_filename, "wb") as f:
                pkl.dump(self.data, f)
            log_debug(f"Saved {self.label} results to cache: {cache_filename}")

    def __enter__(self):
        self.data = self._load_from_cache()
        if self.data is None:
            self.data = {}
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.exists:
            self._save_to_cache()