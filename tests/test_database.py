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
    assert db.getNumberOfDocuments() == 3

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
        "edgelistid": ObjectId("58dbe045ef677d54224a01da"),
        "detectionmethod": "Infomap",
        "r": -1,
        "N": 10,
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


"""
Test the extractDocumentGivenId function
input - an ObjectId
output - a document

Tests:
        - Does it return the right thing given a correct id?
        - Does it fail sensibly given an incorrect id?
        - Does it fail sensibly given a super-incorrect id?

Again, a lightweight wrapper around a mocked method, so can't test too hard.
"""


def test_database_extractdocumentgivenid_normal(mock_database):
    """Assert that the database returns the right document given an id."""
    db = proteinnetworks.database.Database(password="bla")
    doc = db.extractDocumentGivenId("58dbe03fef677d54224a01d9")
    assert doc
