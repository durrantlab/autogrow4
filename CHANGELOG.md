Changes
=======

4.0.1
-----

* Updated `$PATH/autogrow4/accessory_scripts/convert_vina_docked_pdbqt_to_pdbs.py`
  to convert files that label atoms as `HETATM`, Previously, this script only converted
  atoms that were labeled as `ATOM`.
* Updated `$PATH/autogrow4/autogrow/operators/convert_files/conversion_to_3d.py`
  to remove `None` objects which failed to be imported into RDKit. This prevents
  the `None` objects from causing errors in `def convert_single_sdf_to_pdb`.
* Add `raise Exception` in `$PATH/autogrow4/autogrow/operators/convert_files/conversion_to_3d.py`
  with debugging instructions if Gypsum-Dl produced are no SDF files. `raise Exception`
  added to `def convert_sdf_to_pdbs`.
* Added a AutoGrow4 citation to the print statement and `RunAutogrow.py`.
* Added `$PATH/autogrow4/CHANGELOG.MD`.

4.0.0
-----

The original version described in:

Spiegel, J.O., Durrant, J.D. AutoGrow4: an open-source genetic algorithm
for de novo drug design and lead optimization. J Cheminform 12, 25 (2020).
[doi: 10.1186/s13321-020-00429-4]
