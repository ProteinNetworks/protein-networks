"""
Unit tests for the Database class.

Functions:
__init__ x
extractEdgelist
depositEdgelist
extractPDBFile
extractPartition
extractDocumentGivenId
depositPartition o
getNumberOfDocuments x
"""

import proteinnetworks.database
from bson.objectid import ObjectId
import pytest


def test_database_initialise(mock_database):
    """
    Meta-test, to check that pytest and Travis are behaving.

    As the MongoClient interface is mocked for these unit tests,
    the __init__ function doesn't really do anything.
    """
    db = proteinnetworks.database.Database(password="bla")
    assert db


def test_database_getnumberofdocuments(mock_database):
    """
    Test the getNumberOfDocuments function.

    Again, this is a simple wrapper on a mocked class, so hard
    to test.
    """
    db = proteinnetworks.database.Database(password="bla")
    assert db.getNumberOfDocuments() == 4

"""
Test the depositPartition function.

Function inputs: pdbref, edgelistid, detectionmethod, r, N, data.
Function outputs: _id

What is the expected behaviour:
    - When all the arguments are correct? x
    - When the arguments are correct but the partition already
        exists in the database? x
    - When the arguments don't make sense? (Lots of options for this)
"""


def test_database_depositpartition_success(mock_database):
    """Test the depositPartition function when all arguments are correct."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        "pdbref": "1ubq",
        "edgelistid": ObjectId("58dbe03fef677d54224a01da"),
        "detectionmethod": "Infomap",
        "r": -1,
        "N": 5,
        "data": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    }
    resultId = db.depositPartition(**depositionArgs)
    assert type(resultId) == ObjectId


def test_database_depositpartition_doc_already_present(mock_database):
    """Test what happens if you try to add a document that's already there."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': '1ubq',
        'data': [[
            3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4,
            4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 4, 4, 4, 2, 2, 2, 2,
            2, 2, 2, 2
        ], [
            37, 38, 35, 34, 68, 70, 71, 74, 75, 76, 73, 72, 69, 42, 36, 45, 39, 12,
            10, 17, 5, 8, 2, 9, 62, 60, 57, 65, 61, 58, 59, 66, 63, 64, 67, 23, 28,
            26, 29, 21, 19, 22, 47, 49, 46, 55, 56, 54, 53, 50, 14, 11, 18, 7, 6,
            3, 15, 13, 1, 16, 4, 44, 43, 41, 40, 52, 48, 51, 20, 27, 25, 24, 30,
            31, 32, 33
        ]],
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_malformed_pdb(mock_database):
    """Assert that depositPartition handles bad PDB references correctly."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "blablabla",
        'data': [1, 2, 3],
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_edgelistid_invalid(mock_database):
    """Assert that an malformed ObjectId is rejected correctly."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "1ubq",
        'data': [1, 2, 3],
        'N': 10,
        'edgelistid': 'bla',
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_edgelistid_incorrect(mock_database):
    """
    Test valid, but incorrect ObjectIDs.

    Assert that if the ObjectId doesn't correspond to an edgelist
    with the same PDB reference as the partition, the system throws an exception.
    """
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "1ubq",
        'data': [1, 2, 3],
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_detectionmethod_invalid(mock_database):
    """Assert that if "detectionmethod" isn't "AFG" or "Infomap" then an exception is raised."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 2, 3],
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'bla',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_r_invalid(mock_database):
    """Test that invalid r values (aka not a float) are treated properly."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 2, 3],
        'N': -1,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': "ten"
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_N_invalid(mock_database):
    """Test that invalid N values (aka not an integer) are treated properly."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 2, 3],
        'N': "five",
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_r_and_N_not_both_missing(mock_database):
    """Test that if r and N are simultaneously -1 an error is raised."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 2, 3],
        'N': -1,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_r_and_N_not_both_present(mock_database):
    """Test that if r and N are simultaneously not -1 an error is raised."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 2, 3],
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': 1.5
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_data_invalid(mock_database):
    """Test that if the partition isn't an array of ints, an error is raised."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': "bla",
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_data_incorrect(mock_database):
    """
    Test that the system rejects incorrect partitions.

    Partitions should have at least one value in 1,...,M where M is the
    number of unique communities.
    """
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 1, 3, 4],
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


"""
Test the extractDocumentGivenId function
input - an ObjectId
output - a document

Tests:
        - Does it return the right thing given a correct id?
        - Does it fail sensibly given an incorrect id? (i.e one that doesn't exist in the database)
        - Does it fail sensibly given a super-incorrect id?

Again, a lightweight wrapper around a mocked method, so can't test too hard.
"""


def test_database_extractdocumentgivenid_normal(mock_database):
    """Assert that the database returns the right document given an id."""
    db = proteinnetworks.database.Database(password="bla")
    doc = db.extractDocumentGivenId("58dbe03fef677d54224a01d9")
    assert doc


def test_database_extractdocumentgivenid_missing(mock_database):
    """Assert that the database returns None if the document to be found is missing."""
    db = proteinnetworks.database.Database(password="bla")
    doc = db.extractDocumentGivenId("58dbe03fef677d54224a01d7")
    assert not doc


def test_database_extractdocumentgivenid_malformed(mock_database):
    """Assert that the database returns an IOError if the id is malformed."""
    db = proteinnetworks.database.Database(password="bla")
    with pytest.raises(IOError):
        db.extractDocumentGivenId("bla")
