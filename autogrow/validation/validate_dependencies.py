from autogrow.utils.shellcmds import determine_bash_timeout_vs_gtimeout


def validate_dependencies():
    """
    This function will try to import all the installed dependencies that will be
    used in Autogrow. If it fails to import it will raise an ImportError
    """

    # Check Bash Timeout function (There's a difference between MacOS and linux)
    # Linux uses timeout while MacOS uses gtimeout
    timeout_option = determine_bash_timeout_vs_gtimeout()
    if timeout_option not in ["timeout", "gtimeout"]:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
        Autogrow or you may need to execute through Bash."
        )

    try:
        import rdkit
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
        from rdkit.Chem.Draw import PrepareMolForDrawing
        from rdkit.Chem import rdFMCS
        from rdkit.Chem import FilterCatalog
        from rdkit.Chem.FilterCatalog import FilterCatalogParams
        import rdkit.Chem.Lipinski as Lipinski
        import rdkit.Chem.Crippen as Crippen
        import rdkit.Chem.Descriptors as Descriptors
        import rdkit.Chem.MolSurf as MolSurf

    except Exception as e:
        raise ImportError("You need to install rdkit and its dependencies.") from e

    try:
        import numpy
    except Exception as e:
        raise ImportError("You need to install numpy and its dependencies.") from e

    try:
        from scipy.cluster.vq import kmeans2
    except Exception as e:
        raise ImportError("You need to install scipy and its dependencies.") from e
