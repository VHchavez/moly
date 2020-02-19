"""
Unit and regression test for the molely package.
"""

# Import package, test suite, and other packages as needed
import molely
import pytest
import sys

def test_molely_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "molely" in sys.modules
