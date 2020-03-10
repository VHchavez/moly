"""
moly
A short description of the project.
"""

# Add imports here
from .figure.figure import Figure
from .advanced.cubeprop import *
from qcelemental.models.molecule import Molecule

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

