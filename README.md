# AutoGrow4

AutoGrow4 is an evolutionary algorithm that optimizes candidate ligands for
predicted binding affinity and other drug-like properties. Though no
substitute for the medicinal chemist, AutoGrow4 attempts to introduce some
chemical intuition into the automated optimization process.

AutoGrow4 takes advantage of recent advancements in multiprocessing,
cheminformatics, and docking. It handles all ligand manipulations using SMILES
strings via the RDKit API. Docking programs require 3D small-molecule models,
so AutoGrow4 use Gypsum-DL to convert the SMILES to 3D.

AutoGrow4 was designed with modularity in mind. New options allow users to
control drug-likeness, access new reaction libraries, use additional docking
programs, rank AutoGrow4-generated compounds using alternate schemes, and
advance candidate compounds via several selection algorithms. The codebase has
been redesigned so users can easily add their own custom options as well.
Please see the tutorial (`/autogrow4/tutorial/TUTORIAL.md`), which describes
how to expand these various options.

The default fragment libraries included with AutoGrow4 were derived from a
subset of the ZINC database (https://zinc.docking.org/). We thank ZINC for
allowing us to distribute these fragment libraries to AutoGrow4 users.

## Getting Started with AutoGrow4

Use the `RunAutogrow.py` script to run AutoGrow4. The script supports both
command-line parameters and a JSON parameter file. A detailed tutorial
describing AutoGrow4's dependencies and general use is located at
`/autogrow4/tutorial/TUTORIAL.md`. Additional descriptions of all AutoGrow4
parameters can be obtained by running: `python RunAutogrow.py -h`.

We strongly recommend running AutoGrow4 via Docker using
`/autogrow4/docker/autogrow_in_docker.py`. See the tutorial at
`/autogrow4/tutorial/TUTORIAL.md` for more details.

## Dependencies Notes

AutoGrow4 version 4.0.1 has been tested to run using the following dependencies:
  ```python
  >>> rdkit.__version__
  '2019.03.4'
  >>> numpy.__version__
  '1.17.2'
  >>> scipy.__version__
  '1.4.1'
  >>> matplotlib.__version__
  '3.1.1'
  ```

If you are unable to run AutoGrow4, please try running AutoGrow4 in an environment
using these releases these dependencies, or the Docker version of AutoGrow4 which
automatically installs dependencies verified to work with AutoGrow4. If you discover
AutoGrow4 is no longer compatible with current library releases, please contact us 
and we will attempt to correct the code.

## Developer Note

Dimorphite is Version 1.2.3 with the citation print statements silenced.
Please remember to cite Dimorphite-DL: Ropp PJ, Kaminsky JC,
Yablonski S, Durrant JD (2019) Dimorphite-DL: An open-source program for
enumerating the ionization states of drug-like small molecules. J Cheminform
11:14. doi:10.1186/s13321-019-0336-9.

Running AutoGrow4 using Docker uses the most current versions of its dependencies
as of May 2020. These have been tested to be compatible with AutoGrow4. Future
development of AutoGrow4 should test the most current versions of these dependencies
and update `$PATH/docker/Dockerfile` appropriately.