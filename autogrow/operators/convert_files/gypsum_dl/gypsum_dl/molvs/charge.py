# -*- coding: utf-8 -*-
"""
molvs.charge
~~~~~~~~~~~~

This module implements tools for manipulating charges on molecules. In particular, :class:`~molvs.charge.Reionizer`,
which competitively reionizes acids such that the strongest acids ionize first, and :class:`~molvs.charge.Uncharger`,
which attempts to neutralize ionized acids and bases on a molecule.

"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import copy
import logging

from rdkit import Chem

from .utils import memoized_property


log = logging.getLogger(__name__)


class AcidBasePair(object):
    """An acid and its conjugate base, defined by SMARTS.

    A strength-ordered list of AcidBasePairs can be used to ensure the strongest acids in a molecule ionize first.
    """

    def __init__(self, name, acid, base):
        """Initialize an AcidBasePair with the following parameters:

        :param string name: A name for this AcidBasePair.
        :param string acid: SMARTS pattern for the protonated acid.
        :param string base: SMARTS pattern for the conjugate ionized base.
        """
        log.debug("Initializing AcidBasePair: %s", name)
        self.name = name
        self.acid_str = acid
        self.base_str = base

    @memoized_property
    def acid(self):
        log.debug("Loading AcidBasePair acid: %s", self.name)
        return Chem.MolFromSmarts(self.acid_str)

    @memoized_property
    def base(self):
        log.debug("Loading AcidBasePair base: %s", self.name)
        return Chem.MolFromSmarts(self.base_str)

    def __repr__(self):
        return "AcidBasePair({!r}, {!r}, {!r})".format(
            self.name, self.acid_str, self.base_str
        )

    def __str__(self):
        return self.name


#: The default list of AcidBasePairs, sorted from strongest to weakest. This list is derived from the Food and Drug
#: Administration Substance Registration System Standard Operating Procedure guide.
ACID_BASE_PAIRS = (
    AcidBasePair("-OSO3H", "OS(=O)(=O)[OH]", "OS(=O)(=O)[O-]"),
    AcidBasePair("–SO3H", "[!O]S(=O)(=O)[OH]", "[!O]S(=O)(=O)[O-]"),
    AcidBasePair("-OSO2H", "O[SD3](=O)[OH]", "O[SD3](=O)[O-]"),
    AcidBasePair("-SO2H", "[!O][SD3](=O)[OH]", "[!O][SD3](=O)[O-]"),
    AcidBasePair("-OPO3H2", "OP(=O)([OH])[OH]", "OP(=O)([OH])[O-]"),
    AcidBasePair("-PO3H2", "[!O]P(=O)([OH])[OH]", "[!O]P(=O)([OH])[O-]"),
    AcidBasePair("-CO2H", "C(=O)[OH]", "C(=O)[O-]"),
    AcidBasePair("thiophenol", "c[SH]", "c[S-]"),
    AcidBasePair("(-OPO3H)-", "OP(=O)([O-])[OH]", "OP(=O)([O-])[O-]"),
    AcidBasePair("(-PO3H)-", "[!O]P(=O)([O-])[OH]", "[!O]P(=O)([O-])[O-]"),
    AcidBasePair("phthalimide", "O=C2c1ccccc1C(=O)[NH]2", "O=C2c1ccccc1C(=O)[N-]2"),
    AcidBasePair("CO3H (peracetyl)", "C(=O)O[OH]", "C(=O)O[O-]"),
    AcidBasePair("alpha-carbon-hydrogen-nitro group", "O=N(O)[CH]", "O=N(O)[C-]"),
    AcidBasePair("-SO2NH2", "S(=O)(=O)[NH2]", "S(=O)(=O)[NH-]"),
    AcidBasePair("-OBO2H2", "OB([OH])[OH]", "OB([OH])[O-]"),
    AcidBasePair("-BO2H2", "[!O]B([OH])[OH]", "[!O]B([OH])[O-]"),
    AcidBasePair("phenol", "c[OH]", "c[O-]"),
    AcidBasePair("SH (aliphatic)", "C[SH]", "C[S-]"),
    AcidBasePair("(-OBO2H)-", "OB([O-])[OH]", "OB([O-])[O-]"),
    AcidBasePair("(-BO2H)-", "[!O]B([O-])[OH]", "[!O]B([O-])[O-]"),
    AcidBasePair("cyclopentadiene", "C1=CC=C[CH2]1", "c1ccc[cH-]1"),
    AcidBasePair("-CONH2", "C(=O)[NH2]", "C(=O)[NH-]"),
    AcidBasePair("imidazole", "c1cnc[nH]1", "c1cnc[n-]1"),
    AcidBasePair("-OH (aliphatic alcohol)", "[CX4][OH]", "[CX4][O-]"),
    AcidBasePair(
        "alpha-carbon-hydrogen-keto group", "O=C([!O])[C!H0+0]", "O=C([!O])[C-]"
    ),
    AcidBasePair(
        "alpha-carbon-hydrogen-acetyl ester group", "OC(=O)[C!H0+0]", "OC(=O)[C-]"
    ),
    AcidBasePair("sp carbon hydrogen", "C#[CH]", "C#[C-]"),
    AcidBasePair(
        "alpha-carbon-hydrogen-sulfone group", "CS(=O)(=O)[C!H0+0]", "CS(=O)(=O)[C-]"
    ),
    AcidBasePair(
        "alpha-carbon-hydrogen-sulfoxide group", "C[SD3](=O)[C!H0+0]", "C[SD3](=O)[C-]"
    ),
    AcidBasePair("-NH2", "[CX4][NH2]", "[CX4][NH-]"),
    AcidBasePair("benzyl hydrogen", "c[CX4H2]", "c[CX3H-]"),
    AcidBasePair("sp2-carbon hydrogen", "[CX3]=[CX3!H0+0]", "[CX3]=[CX2-]"),
    AcidBasePair("sp3-carbon hydrogen", "[CX4!H0+0]", "[CX3-]"),
)


class ChargeCorrection(object):
    """An atom that should have a certain charge applied, defined by a SMARTS pattern."""

    def __init__(self, name, smarts, charge):
        """Initialize a ChargeCorrection with the following parameters:

        :param string name: A name for this ForcedAtomCharge.
        :param string smarts: SMARTS pattern to match. Charge is applied to the first atom.
        :param int charge: The charge to apply.
        """
        log.debug("Initializing ChargeCorrection: %s", name)
        self.name = name
        self.smarts_str = smarts
        self.charge = charge

    @memoized_property
    def smarts(self):
        log.debug("Loading ChargeCorrection smarts: %s", self.name)
        return Chem.MolFromSmarts(self.smarts_str)

    def __repr__(self):
        return "ChargeCorrection({!r}, {!r}, {!r})".format(
            self.name, self.smarts_str, self.charge
        )

    def __str__(self):
        return self.name


#: The default list of ChargeCorrections.
CHARGE_CORRECTIONS = (
    ChargeCorrection("[Li,Na,K]", "[Li,Na,K;X0+0]", 1),
    ChargeCorrection("[Mg,Ca]", "[Mg,Ca;X0+0]", 2),
    ChargeCorrection("[Cl]", "[Cl;X0+0]", -1),
    # TODO: Extend to other incorrectly charged atoms
)


class Reionizer(object):
    """A class to fix charges and reionize a molecule such that the strongest acids ionize first."""

    def __init__(
        self, acid_base_pairs=ACID_BASE_PAIRS, charge_corrections=CHARGE_CORRECTIONS
    ):
        """Initialize a Reionizer with the following parameter:

        :param acid_base_pairs: A list of :class:`AcidBasePairs <molvs.charge.AcidBasePair>` to reionize, sorted from
                                strongest to weakest.
        :param charge_corrections: A list of :class:`ChargeCorrections <molvs.charge.ChargeCorrection>`.
        """
        log.debug("Initializing Reionizer")
        self.acid_base_pairs = acid_base_pairs
        self.charge_corrections = charge_corrections

    def __call__(self, mol):
        """Calling a Reionizer instance like a function is the same as calling its reionize(mol) method."""
        return self.reionize(mol)

    def reionize(self, mol):
        """Enforce charges on certain atoms, then perform competitive reionization.

        First, charge corrections are applied to ensure, for example, that free metals are correctly ionized. Then, if
        a molecule with multiple acid groups is partially ionized, ensure the strongest acids ionize first.

        The algorithm works as follows:

        - Use SMARTS to find the strongest protonated acid and the weakest ionized acid.
        - If the ionized acid is weaker than the protonated acid, swap proton and repeat.

        :param mol: The molecule to reionize.
        :type mol: rdkit.Chem.rdchem.Mol
        :return: The reionized molecule.
        :rtype: rdkit.Chem.rdchem.Mol
        """
        log.debug("Running Reionizer")

        start_charge = Chem.GetFormalCharge(mol)

        # Apply forced charge corrections
        for cc in self.charge_corrections:
            for match in mol.GetSubstructMatches(cc.smarts):
                atom = mol.GetAtomWithIdx(match[0])
                log.info(
                    "Applying charge correction %s (%s %+d)",
                    cc.name,
                    atom.GetSymbol(),
                    cc.charge,
                )
                atom.SetFormalCharge(cc.charge)

        current_charge = Chem.GetFormalCharge(mol)
        charge_diff = Chem.GetFormalCharge(mol) - start_charge
        # If molecule is now neutral, assume everything is now fixed
        # But otherwise, if charge has become more positive, look for additional protonated acid groups to ionize
        if not current_charge == 0:
            while charge_diff > 0:
                ppos, poccur = self._strongest_protonated(mol)
                if ppos is None:
                    break
                log.info(
                    "Ionizing %s to balance previous charge corrections",
                    self.acid_base_pairs[ppos].name,
                )
                patom = mol.GetAtomWithIdx(poccur[-1])
                patom.SetFormalCharge(patom.GetFormalCharge() - 1)
                if patom.GetNumExplicitHs() > 0:
                    patom.SetNumExplicitHs(patom.GetNumExplicitHs() - 1)
                # else:
                patom.UpdatePropertyCache()
                charge_diff -= 1

        already_moved = set()
        while True:
            ppos, poccur = self._strongest_protonated(mol)
            ipos, ioccur = self._weakest_ionized(mol)
            if ioccur and poccur and ppos < ipos:
                if poccur[-1] == ioccur[-1]:
                    # Bad! H wouldn't be moved, resulting in infinite loop.
                    log.warning("Aborted reionization due to unexpected situation")
                    break

                key = tuple(sorted([poccur[-1], ioccur[-1]]))
                if key in already_moved:
                    log.warning(
                        "Aborting reionization to avoid infinite loop due to it being ambiguous where to put a Hydrogen"
                    )
                    break
                already_moved.add(key)

                log.info(
                    "Moved proton from %s to %s",
                    self.acid_base_pairs[ppos].name,
                    self.acid_base_pairs[ipos].name,
                )

                # Remove hydrogen from strongest protonated
                patom = mol.GetAtomWithIdx(poccur[-1])
                patom.SetFormalCharge(patom.GetFormalCharge() - 1)
                # If no implicit Hs to autoremove, and at least 1 explicit H to remove, reduce explicit count by 1
                if patom.GetNumImplicitHs() == 0 and patom.GetNumExplicitHs() > 0:
                    patom.SetNumExplicitHs(patom.GetNumExplicitHs() - 1)
                    # TODO: Remove any chiral label on patom?
                patom.UpdatePropertyCache()

                # Add hydrogen to weakest ionized
                iatom = mol.GetAtomWithIdx(ioccur[-1])
                iatom.SetFormalCharge(iatom.GetFormalCharge() + 1)
                # Increase explicit H count if no implicit, or aromatic N or P, or non default valence state
                if (
                    iatom.GetNoImplicit()
                    or (
                        (patom.GetAtomicNum() == 7 or patom.GetAtomicNum() == 15)
                        and patom.GetIsAromatic()
                    )
                    or iatom.GetTotalValence()
                    not in list(
                        Chem.GetPeriodicTable().GetValenceList(iatom.GetAtomicNum())
                    )
                ):
                    iatom.SetNumExplicitHs(iatom.GetNumExplicitHs() + 1)
                iatom.UpdatePropertyCache()
            else:
                break

        # TODO: Canonical ionization position if multiple equivalent positions?

        Chem.SanitizeMol(mol)
        return mol

    def _strongest_protonated(self, mol):
        for position, pair in enumerate(self.acid_base_pairs):
            for occurrence in mol.GetSubstructMatches(pair.acid):
                return position, occurrence
        return None, None

    def _weakest_ionized(self, mol):
        for position, pair in enumerate(reversed(self.acid_base_pairs)):
            for occurrence in mol.GetSubstructMatches(pair.base):
                return len(self.acid_base_pairs) - position - 1, occurrence
        return None, None


class Uncharger(object):
    """Class for neutralizing charges in a molecule.

    This class uncharges molecules by adding and/or removing hydrogens. In cases where there is a positive charge that
    is not neutralizable, any corresponding negative charge is also preserved.
    """

    def __init__(self, acid_base_pairs=ACID_BASE_PAIRS):
        log.debug("Initializing Uncharger")
        self.acid_base_pairs = acid_base_pairs
        self.nitro = Chem.MolFromSmarts("[!#8][NX3+](=O)[O-]")

    def __call__(self, mol):
        """Calling an Uncharger instance like a function is the same as calling its uncharge(mol) method."""
        return self.uncharge(mol)

    def uncharge(self, mol):
        """Neutralize molecule by adding/removing hydrogens.

        :param mol: The molecule to uncharge.
        :type mol: rdkit.Chem.rdchem.Mol
        :return: The uncharged molecule.
        :rtype: rdkit.Chem.rdchem.Mol
        """
        log.debug("Running Uncharger")
        mol = copy.deepcopy(mol)

        # Neutralize positive charges
        pos_remainder = 0
        neg_count = 0
        for atom in mol.GetAtoms():
            # Remove hydrogen from positive atoms and reduce formal change until neutral or no more hydrogens
            while atom.GetFormalCharge() > 0 and atom.GetNumExplicitHs() > 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                log.info("Removed positive charge")
            chg = atom.GetFormalCharge()
            if chg > 0:
                # Record number of non-neutralizable positive charges
                pos_remainder += chg
            elif chg < 0:
                # Record total number of negative charges
                neg_count += -chg

        # Choose negative charges to leave in order to balance non-neutralizable positive charges
        neg_skip = self._get_neg_skip(mol, pos_remainder)

        # Neutralize remaining negative charges
        for atom in mol.GetAtoms():
            log.info(atom.GetIdx())
            if atom.GetIdx() in neg_skip:
                continue
            # Make sure to stop when neg_count <= pos_remainder, as it is possible that neg_skip is not large enough
            while atom.GetFormalCharge() < 0 and neg_count > pos_remainder:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                atom.SetFormalCharge(atom.GetFormalCharge() + 1)
                neg_count -= 1
                log.info("Removed negative charge")
        return mol

    def _get_neg_skip(self, mol, pos_count):
        """Get negatively charged atoms to skip (up to pos_count)."""
        neg_skip = set()
        if pos_count:
            # Get negative oxygens in charge-separated nitro groups TODO: Any other special cases to skip?
            for occurrence in mol.GetSubstructMatches(self.nitro):
                neg_skip.add(occurrence[-1])
                if len(neg_skip) >= pos_count:
                    return neg_skip
            # Get strongest ionized acids
            for position, pair in enumerate(self.acid_base_pairs):
                for occurrence in mol.GetSubstructMatches(pair.base):
                    neg_skip.add(occurrence[-1])
                    if len(neg_skip) >= pos_count:
                        return neg_skip
        return neg_skip
