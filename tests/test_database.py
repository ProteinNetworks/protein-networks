"""
Unit tests for the Database class.

Functions:
__init__
extractEdgelist
depositEdgelist
extractPDBFile
extractPartition
extractDocumentGivenId
depositPartition
getNumberOfDocuments
"""


import proteinnetworks.database


def test_initialise(mock_database):
    """meta-test, to check that pytest and Travis are behaving."""
    db = proteinnetworks.database.Database(test=True)
    assert db
