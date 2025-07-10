"""
A fragment filter based on molecular weight.
"""
from typing import Any, List, Tuple
from autogrow.plugins.fragment_filters.fragment_filter_base import FragmentFilterBase
from autogrow.config.argument_vars import ArgumentVars
from autogrow.plugins.registry_base import plugin_managers

class MolecularWeightFilter(FragmentFilterBase):
    """
    A filter for fragments based on their molecular weight.
    This filter is enabled by using the `--MolecularWeightFilter` flag.
    The min and max values can be controlled with `--min_fragment_mol_weight`
    and `--max_fragment_mol_weight`.
    """

    def __init__(self):
        super().__init__()
        self.min_mol_weight = None
        self.max_mol_weight = None

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments specific to this filter.
        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
                group name and a list of argument definitions.
        """
        return (
            "Fragment Filters",
            [
                ArgumentVars(
                    name="MolecularWeightFilter",
                    action="store_true",
                    default=False,
                    help="Enable molecular weight filter for fragments used in mutation."
                ),
                ArgumentVars(
                    name="min_fragment_mol_weight",
                    type=float,
                    default=0.0,
                    help="Minimum molecular weight for fragments. Used with --MolecularWeightFilter.",
                ),
                ArgumentVars(
                    name="max_fragment_mol_weight",
                    type=float,
                    default=500.0,
                    help="Maximum molecular weight for fragments. Used with --MolecularWeightFilter.",
                ),
            ],
        )

    def validate(self, params: dict):
        """
        Validate the parameters for this filter.
        Args:
            params (dict): A dictionary of parameters.
        """
        self.min_mol_weight = params.get("min_fragment_mol_weight")
        self.max_mol_weight = params.get("max_fragment_mol_weight")

    def run_filter(self, mol: Any) -> bool:
        """
        Run the molecular weight filter on a single molecule.
        Args:
            mol: An RDKit molecule object.
        Returns:
            bool: True if the molecule's weight is within the specified range, False otherwise.
        """
        if mol is None:
            return False

        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        mol_weight = chemtoolkit.descriptors_exact_mol_wt(mol)

        return self.min_mol_weight <= mol_weight <= self.max_mol_weight