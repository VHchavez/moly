"""
moly
A short description of the project.
"""

# Add imports here
from .figure.figure import Figure
from qcelemental.models.molecule import Molecule
from .psi4 import psi4
from .data.data_molecules import water, buckyball


# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

