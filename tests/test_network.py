"""
Unit tests for the Network class.

Units to be tested:

__init__.
 generateEdgelist
 extractAtomicData
"""
import proteinnetworks.network
import proteinnetworks.database
from bson.objectid import ObjectId
"""Test the __init__() function in the network.py module.

Inputs: pdbref, edgelisttype, hydrogenstatus, scaling, database
Output: a Network() class with members:
    - scaling
    - edgelisttype
    - hydrogenstatus
    - pdbref
    - database
    - edgelist

Input options:
 - Network params match edgelist in database.
 - Network params correct, but no match
 - Network params mangled (many options, but already tested in test_database)
"""


def test_network_init_edgelist_in_database(mock_database):
    """Test the pathway in which the edgelist required is already in the database."""
    db = proteinnetworks.database.Database(password="bla")
    inputArgs = {
        "scaling": 4.5,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": "2vcr",
        "database": db
    }
    pn = proteinnetworks.network.Network(**inputArgs)
    assert pn.edgelist == [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56],
                           [4, 3, 70], [5, 3, 23]]


def test_network_init_edgelist_not_in_database(mock_database):
    """Test the pathway in which the edgelist must be generated anew."""
    db = proteinnetworks.database.Database(password="bla")
    inputArgs = {
        "scaling": 5.5,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": "1ubq",
        "database": db
    }
    pn = proteinnetworks.network.Network(**inputArgs)
    assert type(pn.edgelistid) == ObjectId


"""
Tests for generateEdgelist

Inputs: pdbref, edgelisttype, hydrogenstatus, scaling.

Output: An edgelist

No tests for malformed input, since checks already done.
Only need to test the case when PDB file not found, and must be fetched.
"""


def test_network_init_edgelist_and_pdb_not_in_database(mock_database):
    """Test the pathway in which the pdb file must be fetched."""
    db = proteinnetworks.database.Database(password="bla")
    inputArgs = {
        "scaling": 4.5,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": "4rvs",
        "database": db
    }
    pn = proteinnetworks.network.Network(**inputArgs)
    assert type(pn.edgelistid) == ObjectId


def test_network_init_edgelist_and_pdb_not_in_database_atomic(mock_database):
    """Test the pathway in which the pdb file must be fetched, with atomic networks."""
    db = proteinnetworks.database.Database(password="bla")
    inputArgs = {
        "scaling": 4.5,
        "edgelisttype": "atomic",
        "hydrogenstatus": "noH",
        "pdbref": "4rvs",
        "database": db
    }
    pn = proteinnetworks.network.Network(**inputArgs)
    assert type(pn.edgelistid) == ObjectId
