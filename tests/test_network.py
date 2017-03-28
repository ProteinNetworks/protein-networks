"""
Unit tests for the Network class.

Units to be tested:

__init__.
    note that this includes Database __init__ (is perhaps a sign that the class design is iffy)

generateEdgelist


def extractAtomicData(self, pdbdata):
"""
import proteinnetworks.network
# import pytest


# This will be run for every function: we replace the actual database with a dict
# @pytest.fixture(autouse=True)
# def mock_database(monkeypatch):
#     monkeypatch.delattr("requests.sessions.Session.request")


def test_extractAtomicData(monkeypatch):
    """Test the extractAtomicData() function in the network.py module.

    Function extractAtomicData reads in pdbdata in the form of a list of strings.
    Outputs positions, elements, residues. Should test:
    """
    pdbdata = ["ATOM bla bla bla"]
    pos, elem, res = proteinnetworks.network.extractAtomicData(pdbdata)