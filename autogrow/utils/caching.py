"""Simple caching utility for storing and retrieving pickled data.

This module provides a context manager-based approach to caching data using
pickle files, with automatic loading and saving of cached results.
"""

from typing import Any, Optional
import os
import pickle as pkl
from autogrow.utils.logging import log_debug, log_info


class CacheManager:
    """A context manager for handling cached data storage and retrieval.

    This class implements a simple caching mechanism using pickle files. It can
    be used as a context manager to automatically load cached data on entry and
    save modified data on exit.

    Args:
        label (str): Identifier used in the cache filename
        cache_dir (str): Directory path where cache files will be stored

    Examples:
        >>> with CacheManager('mydata', '/tmp/cache') as cache:
        ...     if not cache.exists:
        ...         cache.data['new_key'] = compute_expensive_operation()
        ...     result = cache.data['new_key']
    """

    def __init__(self, label: str, cache_dir: str):
        """Initializes a new cache manager instance.

        Args:
            label (str): Identifier used in the cache filename
            cache_dir (str): Directory path where cache files will be stored

        The cache manager is initialized with no data (self.data = None) and
            exists = False until a cache file is actually found during load
            operations.
        """
        self.label = label
        self.cache_dir = cache_dir
        self.exists = False
        self.data: Any = None

    def _get_cache_filename(self) -> str:
        """Constructs the full path to the cache file.

        Returns:
            str: Absolute path to the cache file
        """
        return os.path.join(self.cache_dir, f"{self.label}_results_cache.pkl")

    def _load_from_cache(self) -> Optional[Any]:
        """Attempts to load cached data from disk.

        If the cache file exists, loads and returns its contents. Logs the
        loading operation at INFO level.

        Returns:
            Optional[Any]: The cached data if available, None otherwise
        """
        cache_filename = self._get_cache_filename()
        self.exists = os.path.exists(cache_filename)
        if not self.exists:
            return None
        with open(cache_filename, "rb") as f:
            c = pkl.load(f)
            log_info(
                f"Loaded previous {self.label} results from cache: {cache_filename}"
            )
            return c

    def _save_to_cache(self) -> None:
        """Saves the current data to the cache file if it doesn't already exist.

        Only saves if there is data to save (self.data is not None) and the
        cache doesn't already exist. Logs the save operation at DEBUG level.
        """
        if not self.exists and self.data is not None:
            cache_filename = self._get_cache_filename()
            # TODO: What if dir of cache_filename doesn't exist? Happened once.
            with open(cache_filename, "wb") as f:
                pkl.dump(self.data, f)
            log_debug(f"Saved {self.label} results to cache: {cache_filename}")

    def __enter__(self):
        """Context manager entry point that loads or initializes cache data.

        Returns:
            CacheManager: The cache manager instance with loaded data
        """
        self.data = self._load_from_cache()
        if self.data is None:
            self.data = {}
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit point that saves cache data if necessary.

        Args:
            exc_type: The type of the exception that occurred, if any
            exc_val: The instance of the exception that occurred, if any
            exc_tb: The traceback of the exception that occurred, if any
        """
        if not self.exists:
            self._save_to_cache()
