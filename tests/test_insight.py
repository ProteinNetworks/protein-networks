from bson.objectid import ObjectId
import proteinnetworks.partition
import proteinnetworks.insight

"""
Unit tests for the Insight module.

Units

SuperNetwork
    __init__
    fromPartitionId (classmethod, weird)
    draw
    getIsomorphs
    getWeakIsomorphs

getModifiedJaccard
getZScore
generateNullModel
getMCS
getShannonEntropy
getMutualInfo
getNMI
getConductanceFromPartition
getConductanceFromNodeSubset
getModularityFromPartition
getModularityFromAdjacencyMatrix
"""


"""
Tests for SuperNetwork.__init__()
Inputs: inputPartition, level (int, defaults to None)

Outputs: a SuperNetwork

Options: Multi-chain input partition legit
         Single-chain input partiton legit
         Partition doesn't have given level
         Arguments are invalid

"""


def test_supernetwork_init_all_options_valid_no_level(mock_database):
    """
    Checks that valid inputs gives the correct output.

    If the partition is a valid Partition, and no level is given then a Supernetwork is generated.
    I.e. no supernetwork is in the database, and one is created.
    """
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)  # Tested by partition tests
    # superNetwork = proteinnetworks.insight.SuperNetwork(partition)
    assert partition
