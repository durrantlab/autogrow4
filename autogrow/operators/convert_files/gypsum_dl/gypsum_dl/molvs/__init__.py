# -*- coding: utf-8 -*-
"""
MolVS - Molecule Validation and Standardization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MolVS is a python tool built on top of RDKit that performs validation and standardization of chemical structures.

"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

from .standardize import (
    Standardizer,
    standardize_smiles,
    enumerate_tautomers_smiles,
    canonicalize_tautomer_smiles,
)
from .validate import Validator, validate_smiles
from .errors import MolVSError, StandardizeError, ValidateError


__title__ = "MolVS"
__version__ = "0.1.1"
__author__ = "Matt Swain"
__email__ = "m.swain@me.com"
__license__ = "MIT"
__copyright__ = "Copyright 2019 Matt Swain"


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())
