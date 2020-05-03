Changes
=======

WIP (4.0.1)
-----------

* Updated Gypsum-DL to version 1.1.4.
* Updated Dimorphite-DL to version 1.2.3.
* Revised
  `$PATH/autogrow4/accessory_scripts/convert_vina_docked_pdbqt_to_pdbs.py` to
  convert files that label atoms as `HETATM`, Previously, this script only
  converted atoms that were labeled as `ATOM`.
* Revised
  `$PATH/autogrow4/autogrow/operators/convert_files/conversion_to_3d.py` to
  remove `None` objects which failed to be imported into RDKit. This prevents
  the `None` objects from causing errors in `def convert_single_sdf_to_pdb`.
* Add `raise Exception` in
  `$PATH/autogrow4/autogrow/operators/convert_files/conversion_to_3d.py` with
  debugging instructions if Gypsum-Dl produced are no SDF files. `raise
  Exception` added to `def convert_sdf_to_pdbs`.
* Add `raise Exception` in
  `$PATH/autogrow4/docker/autogrow_in_docker.py` to check for `sudo` and administer
  privileges. 
* Add `--override_sudo_admin_privelledges` variable for 
  `$PATH/autogrow4/docker/autogrow_in_docker.py` skip the check for `sudo` and administer
  privileges for docker-compatible OS that do not have such permissions.
* Fix dependency versions installed in `$PATH/autogrow4/docker/Dockerfile`. This
  ensures that AutoGrow4 continues to run in Docker even if dependency updates
  cause compatibility issues.
* Added a AutoGrow4 citation to the print statement and `$PATH/autogrow4/RunAutogrow.py`.
* Revised `./docker/README.md` to clarify docker use.
* Updated the docker scripts. Bug fixes, added `--rm` to the docker command,
  etc.
* Revised `$PATH/autogrow4/autogrow/user_vars.py` to correct minor typos with OS checks.
* Updated `$PATH/autogrow4/README.md` with minor updates and developer notes.
* Added `Dependencies Notes` to `$PATH/autogrow4/README.md`.
* Added `Developer Notes` to `$PATH/autogrow4/docker/README.md`.
* Added `$PATH/autogrow4/CHANGELOG.md`.

4.0.0
-----

The original version described in:

Spiegel, J.O., Durrant, J.D. AutoGrow4: an open-source genetic algorithm for
de novo drug design and lead optimization. J Cheminform 12, 25 (2020). [doi:
10.1186/s13321-020-00429-4]
