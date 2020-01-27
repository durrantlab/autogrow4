#   AutoGrow 4.0
AutoGrow 4.0 is a complete rewrite of the evolutionary algorithm that optimizes candidate ligands for predicted binding affinity and other druglike properties. Though no substitute for the medicinal chemist, AutoGrow attempts to introduce some chemical intuition into the automated optimization process. By carefully crafting chemically feasible druglike molecules, we hope that AutoGrow 4.0 will help supplement the chemist's efforts.

Autogrow 4.0 takes advantage of the recent advancements in multiprocessing, cheminformatics, and docking. Autogrow 4.0 handles all ligand manipulations using 1D SMILES strings using the API RDKit, which is faster and more maintainable than handling manipulations in 3D using PDB's. As docking softwares require 3D versions of the file, we use Gypsum-DL to convert the 1D SMILES to 3D and Dimorphite-DL to handle the protanation of the ligands using biologically relevant pH settings.

Autogrow 4.0 was designed with modularity in mind. We've expanded new user options for drug-likeliness, reaction library options for Crossover, docking executable softwares, ranking, and selection algorithms. All of these as well as scoring functions have been designed to allow the user to add their own custom options. Please see the tutorial for how and where to expand these options. We hope Autogrow will become a living code base.

A detailed tutorial for using AutoGrow can be found in the directory /autogrow4/tutorial/tutorial.md


The default fragment libraries included with AutoGrow were derived from a subset of the ZINC database (https://zinc.docking.org/). We thank ZINC for allowing us to distribute these fragment libraries to AutoGrow users.

# Getting Started with AutoGrow4:
A detailed tutorial of running AutoGrow4 can be found within the tutorial folder at: /autogrow4/tutorial/tutorial.md
The tutorial details all dependencies, and running instructions. Additional details for AutoGrow4 parameters can be obtained by running:
    python RunAutogrow.py -h

AutoGrow4 is executed through the script RunAutogrow.py and supports commandline and JSON file submission for parameters. 
AutoGrow4 can also be run through Docker, using the script: /autogrow4/Docker/autogrow_in_docker.py
Details for running AutoGrow4 using Docker are provided within the tutorial (/autogrow4/tutorial/tutorial.md)


# Developer Note:
Gypsum-DL is Version 1.1.2 with two sys.flush() commands added to the Parallelizer.py script.
    These help to ensure print statements properly output in large MPI runs.
Dimorphite is Version 1.2.2 with the citation print statement commented out.
    This was done to keep the print logs resonably sized in larger AutoGrow4 runs.
    Please remember to cite Dimorphite-DL:
        Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) 
        Dimorphite-DL: An open-source program for enumerating the ionization states of drug-like small molecules. 
        J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.
