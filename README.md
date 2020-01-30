# AutoGrow4

AutoGrow4 is an evolutionary algorithm that optimizes candidate ligands for
predicted binding affinity and other drug-like properties. Though no
substitute for the medicinal chemist, AutoGrow4 attempts to introduce some
chemical intuition into the automated optimization process.

AutoGrow4 takes advantage of recent advancements in multiprocessing,
cheminformatics, and docking. It handles all ligand manipulations using 1D
SMILES strings via the RDKit API. Docking programs require 3D small-molecule
models, so AutoGrow4 use Gypsum-DL to convert the 1D SMILES to 3D.

AutoGrow4 was designed with modularity in mind. New options allow users to
control drug-likeness, access new reaction libraries, use additional docking
programs, rank AutoGrow4-generated compounds using alternate schemes, and
advance candidate compounds via several selection algorithms. The codebase has
been redesigned so users can easily add their own custom options as well.
Please see the tutorial, which describes how to expand these various options.

The default fragment libraries included with AutoGrow4 were derived from a
subset of the ZINC database (https://zinc.docking.org/). We thank ZINC for
allowing us to distribute these fragment libraries to AutoGrow4 users.

# Getting Started with AutoGrow4

Use the `RunAutogrow.py` script to run AutoGrow4. The script supports both
command-line parameters and a JSON parameter file. A detailed tutorial
describing AutoGrow4's dependencies and general use is located at
`/autogrow4/tutorial/tutorial.md`. Additional descriptions of all AutoGrow4
parameters can be obtained by running: `python RunAutogrow.py -h`.

We strongly recommend running AutoGrow4 via Docker using
`/autogrow4/Docker/autogrow_in_docker.py`. See the tutorial at
`/autogrow4/tutorial/tutorial.md` for more details.

# Developer Note

Gypsum-DL is Version 1.1.2 with two `sys.flush()` commands added to the
`Parallelizer.py` script. These minor changes help ensure that print
statements properly output in large MPI runs.

Dimorphite is Version 1.2.2 with the citation print statement commented out.
We commented out the citation to prevent the AutoGrow4 print logs from growing
too large. Please remember to cite Dimorphite-DL: Ropp PJ, Kaminsky JC,
Yablonski S, Durrant JD (2019) Dimorphite-DL: An open-source program for
enumerating the ionization states of drug-like small molecules. J Cheminform
11:14. doi:10.1186/s13321-019-0336-9.
