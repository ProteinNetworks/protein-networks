"""
Unit tests for the Network class.

Units to be tested:

Network()
    __init__.
    getAdjacencyMatrix
    generateEdgelist
    draw x
extractAtomicData
"""
import proteinnetworks.network
import proteinnetworks.database
from bson.objectid import ObjectId
import numpy as np
"""Test the __init__() function in the network.py module.

Inputs: pdbref, edgelisttype, hydrogenstatus, scaling,
        chainref (optional), database (optional)
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
 - All the above, but with a chainref param.
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


def test_network_init_edgelist_in_database_singlechain(mock_database):
    """Test that if the edgelist is already present, it is successfully generated."""
    pass


def test_network_init_edgelist_not_in_database_singlechain(mock_database):
    """Test when the edgelist must be generated from a multi-chain PDB file."""
    pass


def test_network_init_edgelist_singlechain_chainref_invalid(mock_database):
    """Test when the chainref supplied is invalid (i.e. not a single char)."""
    pass


def test_network_init_edgelist_singlechain_chainref_not_found(mock_database):
    """Test when the chainref supplied is not present (i.e. no ATOM lines contain that chain)."""
    pass


"""
Tests for generateEdgelist

Inputs: pdbref, edgelisttype, hydrogenstatus, scaling, chainref (optional).

Output: An edgelist

No tests for malformed input, since checks already done.
Only need to test the case when PDB file not found, and must be fetched.
Should also test chainref
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


"""
Tests for the Network.draw() method.

Assume that the edgelist is well-formed, i.e. the Network class has valid parameters
(as that should be tested elsewhere)

Inputs: None
Outputs: None. Generates a NetworkX figure.
"""


def test_network_draw_method_creates_plot(mock_database, mock_matplotlib):
    """
    Test that a valid Network instance can create png plots with its draw method.

    Currently all this tests is that the draw() method produces no errors: the actual
    plt.show() call is monkeypatched.
    """
    db = proteinnetworks.database.Database(password="bla")
    inputArgs = {
        "scaling": 4.5,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": "2vcr",
        "database": db
    }
    pn = proteinnetworks.network.Network(**inputArgs)
    pn.draw()


"""
Tests for Network.getAdjacencyMatrix()

Assumes a well-formed Network instance.

Inputs: None.
Member Data Used: Edgelist
Returns: the adjacency Matrix.
"""


def test_network_getadjacencymatrix(mock_database):
    """
    Test that a valid network will return the correct adjacency matrix
    (worked out by hand).
    FIXME we have baked-in assumptions about integer weighting
    """
    db = proteinnetworks.database.Database(password="bla")
    inputArgs = {
        "scaling": 4.5,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": "2vcr",
        "database": db
    }
    pn = proteinnetworks.network.Network(**inputArgs)
    """
    edgelist: [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]]
    expected adjacency matrix:
    0   44  40  0   0
    44  0   56  56  0
    40  56  0   70  23
    0   56  70  0   0
    0   0   23   0  0
    """
    adj = pn.getAdjacencyMatrix()
    expected_adj = np.asarray([[0, 44, 40, 0, 0],
                               [44, 0, 56, 56, 0],
                               [40, 56, 0, 70, 23],
                               [0, 56, 70, 0, 0],
                               [0, 0, 23, 0, 0]], dtype=int)
    assert np.array_equal(adj, expected_adj)
