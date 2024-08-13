__version__ = "4.0.3"


def program_info():
    """
    Get the program version number, etc.

    Returns:
    :returns: str program_output: a string for the print of the program information
    """
    return (
        f"\nAutoGrow Version {__version__}\n"
        + " ================== \n"
        + f"If you use AutoGrow {__version__} in your research, please cite the following reference:\n"
        + "Spiegel, J.O., Durrant, J.D. \n"
        + "AutoGrow4: an open-source genetic algorithm "
        + "for de novo drug design and lead optimization. \n"
        + "J Cheminform 12, 25 (2020). \n"
        + "[doi: 10.1186/s13321-020-00429-4]\n"
        + " ================== \n\n"
    )

