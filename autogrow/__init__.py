"""AutoGrow core module providing version and citation information."""

__version__ = "4.0.3"


def program_info() -> str:
    """Gets program version and citation information for AutoGrow.

    Returns:
        str: A formatted string containing the AutoGrow version number and
            complete citation information for the primary AutoGrow publication.
    """
    return (
        f"\nAutoGrow Version {__version__}\n\n"
        + f"If you use AutoGrow {__version__} in your research, please cite the following reference:\n"
        + "Spiegel, J.O., Durrant, J.D. \n"
        + "AutoGrow4: an open-source genetic algorithm "
        + "for de novo drug design and lead optimization. \n"
        + "J Cheminform 12, 25 (2020). \n"
        + "[doi: 10.1186/s13321-020-00429-4]\n"
    )
