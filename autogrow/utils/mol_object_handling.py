"""Utilities for handling and sanitizing RDKit molecule objects.

This module provides functions for sanitizing RDKit molecules, handling
hydrogens, removing atoms, adjusting nitrogen charges, and managing molecular
fragments. It is adapted from Gypsum-DL
(https://github.com/durrantlab/gypsum_dl).
"""

# Copyright 2018 Jacob D. Durrant

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# NOTE: Adapted from Gypsum-DL. See
# https://github.com/durrantlab/gypsum_dl/blob/main/gypsum_dl/MolObjectHandling.py

##### MolObjectHandling.py
import __future__
from typing import Any, List
from autogrow.plugins.registry_base import plugin_managers


def check_sanitization(mol):
    """Sanitizes an RDKit molecule and fixes common valence errors.

    Attempts to sanitize the molecule using RDKit's sanitization operations. If
    standard sanitization fails, attempts to fix common issues like incorrect
    nitrogen charges before retrying sanitization.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to sanitize

    Returns:
        rdkit.Chem.rdchem.Mol: Sanitized molecule, or None if sanitization fails

    Note:
        Includes a nitrogen-fixing step to correct for common RDKit valence
        errors where nitrogens with 4 bonds have incorrect formal charges.
    """
    if mol is None:
        return None

    chemtoolkit = plugin_managers.ChemToolkit.toolkit

    # easiest nearly everything should get through
    try:
        mol, sanitize_msg = chemtoolkit.sanitize_mol(mol, catch_errors=True,)

    except Exception:
        return None

    if sanitize_msg.name == "SANITIZE_NONE":
        return mol

    # try to fix the nitrogen (common problem that 4 bonded Nitrogens improperly lose their + charges)
    mol = nitrogen_charge_adjustment(mol)

    mol, _ = chemtoolkit.sanitize_mol(mol, catch_errors=True,)
    mol, sanitize_msg = chemtoolkit.sanitize_mol(mol, catch_errors=True,)
    if sanitize_msg.name == "SANITIZE_NONE":
        return mol

    # run a  sanitation Filter 1 more time incase something slipped through
    # ie. if there are any forms of sanition which fail ie. KEKULIZE then return None
    mol, sanitize_msg = chemtoolkit.sanitize_mol(mol, catch_errors=True,)

    if sanitize_msg.name != "SANITIZE_NONE":
        return None

    return mol


def handleHs(mol, protanate_step):
    """
    Control hydrogen atom handling in molecules.

    Sanitizes the molecule and manages explicit/implicit hydrogens based on the
    protanation parameter.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to process
        protanate_step (bool): If True, adds implicit hydrogens. If False,
            removes non-explicit hydrogens

    Returns:
        rdkit.Chem.rdchem.Mol: Processed molecule with handled hydrogens, or
            None if processing fails

    Note:
        Protanation can increase SmilesMerge processing time by up to 10x.
    """
    mol = check_sanitization(mol)
    if mol is None:
        # mol failed Sanitation
        return None

    mol = try_deprotanation(mol)
    if mol is None:
        # mol failed deprotanation
        return None

    if protanate_step is True:
        # PROTANTION IS ON
        mol = try_reprotanation(mol)
        if mol is None:
            # mol failed reprotanation
            return None

    return mol


def try_deprotanation(sanitized_mol):
    """
    Remove non-explicit hydrogens from a sanitized molecule.

    Args:
        sanitized_mol (rdkit.Chem.rdchem.Mol): Sanitized molecule to deprotonate

    Returns:
        rdkit.Chem.rdchem.Mol: Molecule with hydrogens removed and re-sanitized,
            or None if process fails
    """
    try:
        chemtoolkit = plugin_managers.ChemToolkit.toolkit
        mol = chemtoolkit.remove_hs(sanitized_mol, sanitize=False)
    except Exception:
        return None

    return check_sanitization(mol)


def try_reprotanation(sanitized_deprotanated_mol):
    """
    Add implicit hydrogens to a sanitized, deprotonated molecule.

    Args:
        sanitized_deprotanated_mol (rdkit.Chem.rdchem.Mol): Sanitized and
            deprotonated molecule

    Returns:
        rdkit.Chem.rdchem.Mol: Molecule with implicit hydrogens added and
            re-sanitized, or None if process fails
    """
    chemtoolkit = plugin_managers.ChemToolkit.toolkit
    if sanitized_deprotanated_mol is None:
        return None
    try:
        mol = chemtoolkit.add_hs(sanitized_deprotanated_mol)
    except Exception:
        mol = None

    return check_sanitization(mol)


def remove_atoms(mol, list_of_idx_to_remove: List[int]):
    """
    Remove specified atoms from a molecule.

    Uses RDKit's EditableMol class to remove atoms from the molecule based on
    their indices.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to modify
        list_of_idx_to_remove (list): List of atom indices to remove

    Returns:
        rdkit.Chem.rdchem.Mol: Modified molecule with specified atoms removed,
            or None if process fails
    """
    if mol is None:
        return None

    try:
        atoms_to_remove = list_of_idx_to_remove
        atoms_to_remove.sort(reverse=True)
    except Exception:
        return None

    chemtoolkit = plugin_managers.ChemToolkit.toolkit

    try:
        return chemtoolkit.remove_atoms(mol, atoms_to_remove)
    except Exception:
        return None


def nitrogen_charge_adjustment(mol: Any):
    """
    Adjust formal charges on four-bonded nitrogen atoms.

    Corrects for cases where nitrogen atoms have four bonds but lack the
    required positive formal charge. Skips aromatic nitrogens.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to process

    Returns:
        rdkit.Chem.rdchem.Mol: Molecule with corrected nitrogen charges, or None
            if process fails

    Note:
        - Aromatic bonds are treated as 1.5 bond count in RDKit
        - Only non-aromatic nitrogens with exactly 4 bonds are modified
    """
    if mol is None:
        return None
    # makes sure its an rdkit obj

    chemtoolkit = plugin_managers.ChemToolkit.toolkit

    try:
        atoms = chemtoolkit.get_atoms(mol)
    except Exception:
        return None

    for atom in atoms:
        if chemtoolkit.get_atomic_num(atom) == 7:
            bonds = [
                chemtoolkit.get_bond_type_as_double(bond)
                for bond in chemtoolkit.get_bonds(atom)
            ]
            # If aromatic skip as we do not want assume the charge.
            if 1.5 in bonds:
                continue
            # GetBondTypeAsDouble prints out 1 for single, 2.0 for double,
            # 3.0 for triple, 1.5 for AROMATIC but if AROMATIC WE WILL SKIP THIS ATOM
            num_bond_sums = sum(bonds)

            # Check if the octet is filled
            if num_bond_sums == 4.0:
                chemtoolkit.set_formal_charge(atom, +1)
    return mol


def check_for_unassigned_atom(mol):
    """
    Check for presence of unassigned atoms (atomic number 0).

    Identifies if the molecule contains any atoms marked as '*' in SMILES
    notation.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to check

    Returns:
        rdkit.Chem.rdchem.Mol: Input molecule if no unassigned atoms found,
            None otherwise
    """
    if mol is None:
        return None

    chemtoolkit = plugin_managers.ChemToolkit.toolkit

    try:
        atoms = chemtoolkit.get_atoms(mol)
    except Exception:
        return None

    for atom in atoms:
        if chemtoolkit.get_atomic_num(atom) == 0:
            return None
    return mol


def handle_frag_check(mol):
    """
    Process molecules with multiple fragments.

    If the molecule contains multiple fragments, returns the largest fragment
    after checking for unassigned atoms.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to check for fragments

    Returns:
        rdkit.Chem.rdchem.Mol: Largest fragment if molecule is fragmented,
            original molecule if not fragmented, None if processing fails
    """
    if mol is None:
        return None

    chemtoolkit = plugin_managers.ChemToolkit.toolkit

    try:
        frags = chemtoolkit.get_mol_frags(mol, as_mols=True, sanitize_frags=False)
    except Exception:
        return None

    if len(frags) == 1:
        return mol

    frag_info_list = []
    frag_index = 0
    for frag in frags:
        # Check for unassigned breaks ie. a '*'
        frag = check_for_unassigned_atom(frag)
        if frag is None:
            frag_index = frag_index + 1
            continue

        num_atoms = chemtoolkit.get_num_atoms(frag)
        frag_info = [frag_index, num_atoms]
        frag_info_list.append(frag_info)
        frag_index = frag_index + 1
    if len(frag_info_list) == 0:
        return None

    # Get the largest Fragment
    frag_info_list.sort(key=lambda x: float(x[-1]), reverse=True)
    largest_frag_idx = frag_info_list[0][0]
    largest_frag = frags[largest_frag_idx]
    return largest_frag
