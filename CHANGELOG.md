Changes
=======

WIP
---

* Improved error handling. (Actually, instead of this, just use latest version of gypsum. Need to bring that over. That should be the update.)


4.0.3
-----

* Updated Gypsum-DL to version 1.1.7.
* Revised `$PATH/autogrow4/autogrow/user_vars.py` to correct for error in
  handling the `--conversion_choice`. Obabel can now be run as the
  `--conversion_choice` without providing an MGLTools path
  (`mgltools_directory`).
* Revised `$PATH/autogrow4/autogrow/user_vars.py` to enforce
  `MGLToolsConversion` as the only conversion option (`--conversion_choice`)
  when scoring with NN1/NN2 (`--scoring_choice`).
* Minor code clean-up to `$PATH/autogrow4/autogrow/user_vars.py`.
* Fixed [user-reported
  bug](https://durrantlab.pitt.edu/forums/topic/autogrow4-bug-report-in-accessory_scripts-fragmenter_of_smi_mol-py/)
  in `fragmenter_of_smi_mol.py`.

4.0.2
-----

* Updated Gypsum-DL to version 1.1.5.
* Updated Dimorphite-DL to version 1.2.4.
* Fixed failure in `$PATH/autogrow4/accessory_scripts/make_lineage_figures.py`
  to detect source compounds when run had `--use_docked_source_compounds` set
  to False. This patched required the added user variable
  `--use_docked_source_compounds`. To maintain back-compatibility the default
  setting is to auto-detect from the `vars.json` file.
* Added optional variable `--purge_previous_pickled_files` to
  `$PATH/autogrow4/accessory_scripts/make_lineage_figures.py` which automates
  removing all pickled files generated by the script. Useful for cleaning up
  files when done processing lineages and debugging. Default is False so it is
  backwards compatible.
* Revised `$PATH/autogrow4/autogrow/user_vars.py` to correct minor typos with
  running custom filters.
* Removed blank line in `$PATH/autogrow4/source_compounds/naphthalene_smiles.smi`.

4.0.1
-----

* Updated Gypsum-DL to version 1.1.4.
* Updated Dimorphite-DL to version 1.2.3.
* Revised
  `$PATH/autogrow4/accessory_scripts/convert_vina_docked_pdbqt_to_pdbs.py` to
  convert files that label atoms as `HETATM`. Previously, this script only
  converted atoms that were labeled as `ATOM`.
* Revised
  `$PATH/autogrow4/autogrow/operators/convert_files/conversion_to_3d.py` to
  remove `None` objects that could not be imported into RDKit. This prevents
  the `None` objects from causing errors in `def convert_single_sdf_to_pdb`.
* Add `raise Exception` in
  `$PATH/autogrow4/autogrow/operators/convert_files/conversion_to_3d.py` with
  debugging instructions if Gypsum-DL produces no SDF files. `raise Exception`
  added to `def convert_sdf_to_pdbs`.
* Add `raise Exception` in `$PATH/autogrow4/docker/autogrow_in_docker.py` to
  check for `sudo` and administrator privileges.
* Add `--override_sudo_admin_privileges` parameter to
  `$PATH/autogrow4/docker/autogrow_in_docker.py`. This parameter skips the
  check for `sudo` and administrator privileges for docker-compatible OS that
  do not have such permissions.
* Fix dependency versions installed in `$PATH/autogrow4/docker/Dockerfile`.
  This ensures that AutoGrow4 continues to run in Docker even if future
  dependency updates cause compatibility issues.
* Added an AutoGrow4 citation to the print statement in
  `$PATH/autogrow4/run_autogrow.py`.
* Revised `./docker/README.md` to clarify docker use.
* Updated the docker scripts. Bug fixes, added `--rm` to the docker command,
  etc.
* Revised `$PATH/autogrow4/autogrow/user_vars.py` to correct minor typos with
  OS checks.
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
