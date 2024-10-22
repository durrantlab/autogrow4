"""Dependency validation utilities for AutoGrow.

This module verifies that all required Python packages are installed and
importable before AutoGrow execution begins. This includes comprehensive checks
for RDKit and its submodules, numpy, and scipy components.
"""


def validate_dependencies() -> None:
    """Validates that all required dependencies are installed.
   
    Attempts to import all dependencies used by AutoGrow including:
    
    RDKit packages:
        - rdkit core
        - Chem and AllChem
        - rdDepictor
        - rdMolDraw2D and PrepareMolForDrawing 
        - rdFMCS
        - FilterCatalog and FilterCatalogParams
        - Lipinski
        - Crippen
        - Descriptors
        - MolSurf
   
    Other required packages:
        - numpy
        - scipy.cluster.vq.kmeans2
   
    Raises:
        ImportError: If any required package fails to import. Message specifies
            which package is missing.
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
