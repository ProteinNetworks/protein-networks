import proteinnetworks.partition
from bson.objectid import ObjectId
"""
Unit tests for the partition module.

Units:
Partition:
    __init__
    generatePartition
treeFileToNestedLists
"""

"""
Tests for Partition.__init__()

inputs: pdbref, edgelistid, detectionmethod, r, N.
outputs: a Partition

Options: - params legit, partition in database
         - params malformed (will be picked up by DB class)
         - params legit, partition not in database
"""


def test_partition_init_partition_in_database(mock_database):
    """Test that if the partition is in the DB, it is extracted successfully."""
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    assert type(partition.partitionid) == ObjectId


def test_partition_init_partition_not_in_database(mock_database, mock_subprocess):
    """
    Test that if the partition is not in the DB, it is generated successfully.

    Uses a mocked community detection method.
    """
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 1,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    assert type(partition.partitionid) == ObjectId


"""
Tests for the treeFileToNestedLists function (not a method)

Takes in a path to a treefile, outputs a list of lists.

test: if the treefile doesn't exist we throw a FileNotFoundError
      if the treefile format is unexpected we throw an IOError
      if the treefile is legit we return a list of lists.
"""
