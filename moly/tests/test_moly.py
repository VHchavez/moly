"""
Unit and regression test for the moly package.
"""

# Import package, test suite, and other packages as needed
import moly
import pytest
import sys

def test_moly_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "moly" in sys.modules
