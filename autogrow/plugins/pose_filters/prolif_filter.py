"""
ProLIF Filter for filtering following interaction fingerprints.

This module uses the ProLIF library to get interaction fingerprints. To pass the ProLIF filter, a molecule must have
at least one interaction with a given receptor.

ProLIF library is available at https://github.com/chemosim-lab/ProLIF
"""
import __future__

from autogrow.plugins.pose_filters import PoseFilterBase
import prolif
import rdkit  # type: ignore
import rdkit.Chem as Chem  # type: ignore
import rdkit.Chem.Descriptors as Descriptors  # type: ignore
from typing import List, Tuple
from autogrow.config.argparser import ArgumentVars

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


class ProLIFFilter(PoseFilterBase):
    """
    ProLIF Filter for filtering following interaction fingerprints.

    This module uses the ProLIF library to get interaction fingerprints. To pass the ProLIF filter, a molecule must have
    at least one interaction with a given receptor.

    ProLIF library is available at https://github.com/chemosim-lab/ProLIF
    """

    def run_filter(self, receptor, docked_cmpd) -> bool:
        """
        Run the ProLIF filter on a given molecule against given receptor.

        This method calculates interaction fingerprints using the ProLIF library against a given receptor to determine
         if a docked molecule has at least interaction regarding a given receptor. It checks Hydrophobic, HBDonor,
         HBAcceptor, PiStacking, Anionic, Cationic, CationPi, and PiCation interactions.

        Args:
            receptor: A rdkit object representing the receptor
            docked_cmpd: A rdkit object representing the docked molecule

        Returns:
            bool: True if the molecule passes the filter (allows up to one
                violation), False otherwise.
        """
        fp = prolif.Fingerprint(
            ["Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking", "Anionic", "Cationic", "CationPi", "PiCation"])

        try:
            prot = prolif.Molecule.from_rdkit(receptor)
            lig = prolif.Molecule.from_rdkit(docked_cmpd)
            ifp = fp.generate(lig, prot, metadata=True)
            df = prolif.to_dataframe({0: ifp}, fp.interactions)

            return df.size > 0
        except Exception as e:
            return False

    def add_arguments(self) -> Tuple[str, List[ArgumentVars]]:
        """
        Add command-line arguments required by the plugin.

        This method defines the command-line arguments specific to the ProLIF Filter.
        It allows users to enable the filter via command-line options.

        Returns:
            Tuple[str, List[ArgumentVars]]: A tuple containing the argument
                group name and a list of ArgumentVars objects defining the
                command-line arguments.
        """
        return (
            "Pose Filters",
            [
                ArgumentVars(
                    name=self.name,
                    action="store_true",
                    default=False,
                    help="ProLIF interaction fingerprints.",
                )
            ],
        )
