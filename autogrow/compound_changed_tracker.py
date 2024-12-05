from dataclasses import dataclass
from typing import List, Optional
from contextlib import contextmanager
import copy


@dataclass
class CompoundState:
    """Captures the state of a Compound's tracked attributes"""

    sdf_path: Optional[str]
    history: List[str]


class CompoundChangeTracker:
    """Context manager to track changes in Compound attributes."""

    def __init__(self, compound, track_history=True, track_sdf_path=True):
        """Initialize with a Compound instance to track."""
        self.compound = compound
        self.initial_state = None
        self.track_history = track_history
        self.track_sdf_path = track_sdf_path

    def _capture_state(self) -> CompoundState:
        """Capture current state of tracked attributes."""
        return CompoundState(
            sdf_path=self.compound.sdf_path,
            history=copy.deepcopy(self.compound.history),
        )

    def _verify_changes(self):
        """Verify that tracked attributes have changed."""
        current_state = self._capture_state()

        assert self.initial_state is not None, "Initial state was not captured."

        if (
            self.track_sdf_path
            and current_state.sdf_path == self.initial_state.sdf_path
        ):
            raise ValueError("sdf_path was not modified in the compound operation.")

        if self.track_history and current_state.history == self.initial_state.history:
            raise ValueError("history was not modified in the compound operation.")

    def __enter__(self):
        """Enter the context and capture initial state."""
        self.initial_state = self._capture_state()
        return self.compound

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the context and verify changes occurred."""
        if exc_type is None:  # Only verify if no exception occurred
            self._verify_changes()


# Helper function to make usage cleaner
@contextmanager
def track_compound_changes(compound, track_history=True, track_sdf_path=True):
    """Track changes to a compound's sdf_path and history attributes.

    Args:
        compound (Compound): The Compound instance to track.
        track_history (bool): Whether to track changes to the history attribute.
        track_sdf_path (bool): Whether to track changes to the sdf_path attribute.
    
    Usage:
        with track_compound_changes(my_compound):
            # Perform operations that should modify sdf_path or history
            my_compound.sdf_path = "new_path.sdf"
            my_compound.history.append("Operation performed")
    """
    tracker = CompoundChangeTracker(compound, track_history, track_sdf_path)
    with tracker:
        yield compound
