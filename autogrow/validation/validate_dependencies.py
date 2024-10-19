def validate_dependencies() -> None:
    """
    This function will try to import all the installed dependencies that will be
    used in Autogrow. If it fails to import it will raise an ImportError
    """

    try:
        import rdkit  # type: ignore
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem  # type: ignore
        from rdkit.Chem import rdDepictor  # type: ignore
        from rdkit.Chem.Draw import rdMolDraw2D  # type: ignore
        from rdkit.Chem.Draw import PrepareMolForDrawing  # type: ignore
        from rdkit.Chem import rdFMCS  # type: ignore
        from rdkit.Chem import FilterCatalog  # type: ignore
        from rdkit.Chem.FilterCatalog import FilterCatalogParams  # type: ignore
        import rdkit.Chem.Lipinski as Lipinski  # type: ignore
        import rdkit.Chem.Crippen as Crippen  # type: ignore
        import rdkit.Chem.Descriptors as Descriptors  # type: ignore
        import rdkit.Chem.MolSurf as MolSurf  # type: ignore

    except Exception as e:
        raise ImportError("You need to install rdkit and its dependencies.") from e

    try:
        import numpy
    except Exception as e:
        raise ImportError("You need to install numpy and its dependencies.") from e

    try:
        from scipy.cluster.vq import kmeans2  # type: ignore
    except Exception as e:
        raise ImportError("You need to install scipy and its dependencies.") from e
