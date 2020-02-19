"""
moly
A short description of the project.
"""

# Add imports here
from .molecule.molecule_factory import molecule_factory
from .figure.figure import Figure


# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions


