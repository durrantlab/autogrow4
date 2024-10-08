from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.mutation import MutationBase


class FragmentAddition(MutationBase):
    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """Add command-line arguments required by the plugin."""
        # TODO: Update
        return (
            "Selectors",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="Select compounds to advance to the next generation per a ranked selector. The Rank option is a non-redundant selector. Do not use Rank_Selector for small runs as there is potential that the number of desired ligands exceed the number of ligands to chose from.",  # TODO: Add more detail here.
                )
            ],
        )

    def validate(self, params: dict):
        """Validate the provided arguments."""
        pass
