# AutoGrow5

AutoGrow5 is an evolutionary algorithm that optimizes candidate ligands for
predicted binding affinity and other drug-like properties. While it cannot
replace the expertise of medicinal chemists, AutoGrow5 does seek to bring
chemical intuition to the automated optimization process.

## Overview

AutoGrow5 uses an evolutionary algorithm to generate novel compounds with
desirable properties. Starting from a population of seed molecules, it
iteratively applies genetic operators to create new generations of compounds:

- **Mutation**: Modifies molecules (e.g., by adding fragments)
- **Crossover**: Combines parts of two parent molecules to create offspring
- **Elitism**: Advances the best-performing compounds from each generation to
  the next

This evolutionary process guides optimization toward molecules with better
predicted binding affinity and drug-like characteristics.

## Modular Architecture

AutoGrow5 was designed with modularity as a core principle. Its plugin-based
architecture allows extensive customization of the drug discovery workflow:

- **Drug-likeness filtering**: Control compound quality through various
  SMILES-based filters (Lipinski, Ghose, PAINS)
- **Reaction libraries**: Access different libraries for mutation and crossover
  operations
- **Docking programs**: Choose from various options (AutoDock Vina, GNINA)
- **Selection algorithms**: Use different compound selection methods (Rank,
  Tournament, Roulette)
- **Pose filtering**: Apply ProLIF-based interaction filters and other methods
- **Rescoring**: Use ligand efficiency and other rescoring approaches

The codebase has also been redesigned to make it easy for users to add custom
options.

## Fragment Libraries

The default fragment libraries included with AutoGrow5 (used by the default
fragment addition mutation operator) are derived from two sources: the NCI/DTP
Open Chemicals Repository and the original AutoGrow4 libraries (a subset of the
ZINC database). This initial compound set was subjected to rigorous filtering to
ensure drug-likeness, including constraints on molecular weight, rotatable
bonds, and cLogP. In some cases, the libraries were further augmented using in
silico synthesis, retrosynthesis, and stereoisomer variation to ensure a diverse
set of small, extensible fragments. See
`autogrow/plugins/mutation/reaction_libraries/all_rxns/AUTOGROW5_MOIETY_LIB.md`
for full details.

We thank ZINC (https://zinc.docking.org/) for allowing us to distribute these
fragment libraries to AutoGrow5 users.

## Citation

If you use AutoGrow version 5.0 in your research, please cite:

Spiegel, J.O., Durrant, J.D. AutoGrow4: an open-source genetic algorithm for
de novo drug design and lead optimization. J Cheminform 12, 25 (2020). doi:
10.1186/s13321-020-00429-4

## Getting Started

Use the `run_autogrow.py` script to run AutoGrow5. The script supports both
command-line parameters and JSON parameter files. To see all available plugins
and parameters, run:

```bash
python run_autogrow.py -h
```

## Dependencies

Python dependencies for AutoGrow5 are specified in the `environment.yml` file.
AutoGrow5 also requires some third-party software such as OpenBabel for file
conversions. Some accessory scripts for plotting and analysis may require
`matplotlib`, `pandas`, `numpy`, `scikit-learn`, etc.

We have tested AutoGrow5 on macOS and Linux. We expect it can also run run on
Windows, though Windows is not officially supported. If you discover
compatibility issues with current library releases, please contact us, and we
will work to correct the code.

## Environment Setup

To create a Python environment using conda:

1. `conda env create -f environment.yml`
2. `conda activate AutoGrow5`
