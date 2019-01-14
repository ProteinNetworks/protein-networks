from bson.objectid import ObjectId
import proteinnetworks
import pytest
import numpy as np
import networkx as nx
import math

from scipy.linalg import block_diag

"""
Unit tests for the Insight module.

Units

SuperNetwork
    __init__
    fromPartitionId (classmethod, weird)
    draw
    getIsomorphs
    getWeakIsomorphs

SuperNetworkNullModel

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
edgelistToGraph
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
    superNetwork = proteinnetworks.insight.SuperNetwork(partition)
    # Need to check that the partition is as expected
    assert superNetwork


def test_supernetwork_init_all_options_valid_with_level(mock_database):
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
    superNetwork = proteinnetworks.insight.SuperNetwork(partition, level=0)
    # Need to check that the partition is as expected
    assert superNetwork


def test_supernetwork_init_no_level_no_PFAM(mock_database):
    """
    Test that if no partition level is provided, and no PFAM entry can be found, then
    an error is thrown.
    """
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '2vcr',
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)  # Tested by partition tests
    with pytest.raises(ValueError):
        superNetwork = proteinnetworks.insight.SuperNetwork(partition)


def test_supernetwork_init_malformed_level(mock_database):
    """
    Test that if a non-int level is passed to the supernetwork then a TypeError is thrown. 
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
    with pytest.raises(TypeError):
        superNetwork = proteinnetworks.insight.SuperNetwork(partition, level="bla")
    

def test_supernetwork_init_invalid_level(mock_database):
    """
    Test that if a level which doesn't correspond to the partition level is passed then
    an IndexError is thrown.
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
    with pytest.raises(IndexError):
        superNetwork = proteinnetworks.insight.SuperNetwork(partition, level=-10)
    # assert superNetwork
    with pytest.raises(IndexError):
        superNetwork = proteinnetworks.insight.SuperNetwork(partition, level=10)


# def test_supernetwork_retrieved_from_database(mock_database):
#     """
#     Test that if a supernetwork has been stored in the database it can be retrieved correctly.
#     """
#     db = proteinnetworks.database.Database(password="bla")
#     partitionArgs = {
#         'pdbref': '2vcr',
#         'N': 10,
#         'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
#         'detectionmethod': 'Infomap',
#         'r': -1,
#         "database": db
#     }
#     partition = proteinnetworks.partition.Partition(**partitionArgs)  # Tested by partition tests
#     # with pytest.raises(ValueError):
#     superNetwork = proteinnetworks.insight.SuperNetwork(partition, level=0)
#     assert superNetwork

"""
tests for getModifiedJaccard()

inputs:
An expected Array and a generated Array (both numpy arrays)

returns:
the modified jaccard between the two arrays.

tests:
expectedArray and/or generatedArray not numpy arrays (e.g. lists or strings)
expectedArray and/or generatedArray of length 0

the case with overlap of one domain
the case with some overlap of multiple domains
the case with total overlap
multiple expected domains
"""

def test_getmodifiedjaccard_invalidinputs_wrongtype_list():
    """
    Test that non-numpy arrays are rejected
    """
    expected = [1]*40 + [2]*20 + [1]*40
    generated = [1]*40 + [2]*20 + [1]*40
    with pytest.raises(TypeError):
        jaccard = proteinnetworks.insight.getModifiedJaccard(expected, generated)


def test_getmodifiedjaccard_invalidinputs_wrongtype_string():
    """
    Second test that numpy arrays of non-ints (in this case, strings) are rejected
    """
    expected = np.asarray(["a"]*40 + ["b"]*20 + ["a"]*40)
    generated = np.asarray(["a"]*40 + ["b"]*20 + ["c"]*40)
    with pytest.raises(TypeError):
        jaccard = proteinnetworks.insight.getModifiedJaccard(expected, generated)


def test_getmodifiedjaccard_invalidinputs_wrongtype_floats():
    """
    Second test that numpy arrays of non-ints (in this case, floats) are rejected
    """
    expected = np.asarray([0.5]*40 + [1]*20 + [0.5]*40)
    generated = np.asarray([0.5]*40 + [1.0]*20 + [1.5]*40)
    with pytest.raises(TypeError):
        jaccard = proteinnetworks.insight.getModifiedJaccard(expected, generated)


def test_getmodifiedjaccard_invalidinputs_emptyarrays():
    """
    Test that zero-length numpy arrays are rejected
    """
    expected = np.asarray([])
    generated = np.asarray([])
    with pytest.raises(ValueError):
        jaccard = proteinnetworks.insight.getModifiedJaccard(expected, generated)


def test_getmodifiedjaccard_invalidinputs_differentarrays():
    """
    Test that different-length inputs are rejected
    """
    expected = np.asarray([1]*40 + [2]*20 + [1]*40)
    generated = np.asarray([1]*40 + [2]*20)
    with pytest.raises(ValueError):
        jaccard = proteinnetworks.insight.getModifiedJaccard(expected, generated)


def test_getmodifiedjaccard_perfectoverlap():
    """
    Check that if the PFAM array has perfect overlap with one domain, a score of 1.0 is returned. 
    """
    expected = np.asarray([1]*40 + [2]*20 + [1]*40)
    generated = np.asarray([1]*40 + [2]*20 + [3]*20 + [4]*20)
    jaccard = proteinnetworks.insight.getModifiedJaccard(expected, generated)
    assert jaccard == 1


def test_getmodifiedjaccard_partialoverlap_single_domains():
    """
    Check that if a PFAM array is a subset of a community, a score of <1.0 is returned.
    """
    expected = np.asarray([1]*40 + [2]*20 + [1]*40)
    generated = np.asarray([1]*40 + [2]*30 + [3]*10 + [4]*20)
    jaccard = proteinnetworks.insight.getModifiedJaccard(expected, generated)
    # intersection of size 20, union of size 30
    assert jaccard == 2/3


def test_getmodifiedjaccard_partialoverlap_multiple_domains():
    """
    Check that if a PFAM array is a subset of a community, a score of <1.0 is returned.
    """
    expected = np.asarray([1]*40 + [2]*20 + [1]*40)
    generated = np.asarray([1]*40 + [2]*10 + [3]*10 + [4]*40)
    jaccard = proteinnetworks.insight.getModifiedJaccard(expected, generated)
    # mean overlap is 0.5
    assert jaccard == 1/2


def test_getmodifiedjaccard_nopfamdomain():
    """
    Check that if the PFAM array is all 1s (i.e. no PFAM domains) an error is thrown.
    """
    expected = np.asarray([1]*40)
    generated = np.asarray([1]*20 + [2]*20)
    with pytest.raises(ValueError):
        jaccard = proteinnetworks.insight.getModifiedJaccard(expected, generated)


"""
Tests for generateNullModel

given a test partition, generate a partition with the same number of inter-community boundaries
and the same total number of communities

returns a numpy array (hopefully)

tests:
expectedArray not a numpy array (e.g. lists or strings)
expectedArray of length 0
the case with unit-length communities
the case with normal communities

check return is a numpy array.
"""


def test_generatenullmodel_invalidinput_list():
    """Check that lists are rejected."""
    expected = [1]*40
    with pytest.raises(TypeError):
        array = proteinnetworks.insight.generateNullModel(expected)


def test_generatenullmodel_invalidinput_string():
    """Check that arrays of strings are rejected."""
    expected = np.asarray(["bla"]*40)
    with pytest.raises(TypeError):
        array = proteinnetworks.insight.generateNullModel(expected)


def test_generatenullmodel_invalidinput_float():
    """Check that arrays of floats are  rejected."""

    expected = np.asarray([1.5]*40)
    with pytest.raises(TypeError):
        array = proteinnetworks.insight.generateNullModel(expected)


def test_generatenullmodel_invalidinput_zerolengtharray():
    """Check that zero-length arrays are rejected."""
    expected = np.asarray([])
    with pytest.raises(ValueError):
        array = proteinnetworks.insight.generateNullModel(expected)


def test_generatenullmodel_allcommunitiessizeone():
    expected = np.asarray(range(1,41))
    array = proteinnetworks.insight.generateNullModel(expected)
    assert len(set(array)) == 40
    numBoundaries = -1
    prevI = -1
    for i in array:
        if i != prevI:
            numBoundaries += 1
        prevI = i
    assert numBoundaries == 39


def test_generatenullmodel_singlecommunity():
    expected = np.asarray([1]*40)
    array = proteinnetworks.insight.generateNullModel(expected)
    assert len(set(array)) == 1
    numBoundaries = -1
    prevI = -1
    for i in array:
        if i != prevI:
            numBoundaries += 1
        prevI = i
    assert numBoundaries == 0


def test_generatenullmodel_normalcommunities():
    expected = np.asarray([1]*20 + [2]*20 + [3]*20 )
    array = proteinnetworks.insight.generateNullModel(expected)
    assert len(set(array)) == 3
    numBoundaries = -1
    prevI = -1
    for i in array:
        if i != prevI:
            numBoundaries += 1
        prevI = i
    assert numBoundaries == 2


def test_generatenullmodel_checkreturntype():
    expected = np.asarray([1]*40)
    expected = np.asarray([1]*20 + [2]*20 + [3]*20 )
    array = proteinnetworks.insight.generateNullModel(expected)
    assert type(array) == np.ndarray


def test_generatenullmodel_invalidinput_communities_start_from_one():
    """Check that the communities are labelled from 1 to N, without gaps"""
    expected = np.asarray([0]*5+[2]*5)
    with pytest.raises(ValueError):
        array = proteinnetworks.insight.generateNullModel(expected)


"""
tests for getMCS.

Inputs -> two networkx graphs

Outputs -> a smaller graph


Test:
- Invalid input
- Isomorphic graphs
- Graphs of the same size
- Graphs of different sizes
- Graphs with only one node
- Disconnected graphs
"""


def test_getmcs_invalid_input():
    """Test that inputs that aren't networkx undirected graphs are rejected.""" 
    input1 = [[1,2], [2,3], [1,3]]
    input2 = "bla"
    input1 = nx.DiGraph(input1)
    with pytest.raises(TypeError):
        mcs = proteinnetworks.insight.getMCS(input1, input2)


def test_getmcs_isomorphic_graphs():
    input1 = [[1,2], [2,3], [1,3]] # triangle
    input2 = [[2,3], [3,4], [4,2]] # triangle, differently labelled

    input1 = nx.Graph(input1)
    input2 = nx.Graph(input2)
    mcs = proteinnetworks.insight.getMCS(input1, input2)
    assert nx.is_isomorphic(mcs, input1) and nx.is_isomorphic(mcs, input2)


def test_getmcs_perfect_subgraph():

    input1 = [[1,2], [2,3], [1,3]] # triangle
    input2 = [[2,3], [3,4], [4,2], [2,5]] # triangle with an extra node

    input1 = nx.Graph(input1)
    input2 = nx.Graph(input2)
    mcs = proteinnetworks.insight.getMCS(input1, input2)
    assert nx.is_isomorphic(mcs, input1) and not nx.is_isomorphic(mcs, input2)


def test_getmcs_single_node_graph():
    input1 = [[1,2], [2,3], [1,3]] # triangle

    input1 = nx.Graph(input1)
    input2 = nx.Graph()
    input2.add_node(1)
    mcs = proteinnetworks.insight.getMCS(input1, input2)
    assert nx.is_isomorphic(mcs, input2) and not nx.is_isomorphic(mcs, input1)
    assert mcs.number_of_nodes() == 1


def test_getmcs_disconnected_graphs():
    input1 = [[1,2], [2,3], [1,3]] # triangle
    input2 = [[2,3], [3,4], [4,2], [5,6], [6,7], [7,5]] # two triangles
    input1 = nx.Graph(input1)
    input2 = nx.Graph(input2)
    mcs = proteinnetworks.insight.getMCS(input1, input2)
    assert nx.is_isomorphic(mcs, input1) and not nx.is_isomorphic(mcs, input2)


"""
Tests for getZscore

inputs -> two arrays & a number of trials
Output -> a float giving significance

We expect single-community and unit-community partitions to have significance zero (I think)

Tests:

Bad input (negative numTrials)

Single community comparison
Unit community comparison
Normal case comparison

FIXME any what circumstances will this be negative?

"""


def test_getzscore_invalid_number_of_trials():
    """Test that a negative number of trials is rejected."""
    expected = np.asarray([1]*40 + [2]*20 + [1]*40)
    generated = np.asarray([1]*40 + [2]*10 + [3]*10 + [4]*40)
    with pytest.raises(ValueError):
        score = proteinnetworks.insight.getZScore(expected, generated, numTrials=-100)


def test_getzscore_single_community():
    expected = np.asarray([1]*40 + [2]*20 + [1]*40)
    generated = np.asarray([1]*100)
    score = proteinnetworks.insight.getZScore(expected, generated)
    assert score == 0


def test_getzscore_unit_community():
    expected = np.asarray([1]*40 + [2]*20 + [1]*40)
    generated = np.asarray(range(1,101))
    score = proteinnetworks.insight.getZScore(expected, generated)
    assert score == 0


def test_getzscore_normal_communities():
    expected = np.asarray([1]*40 + [2]*20 + [1]*40)
    generated = np.asarray([1]*40 + [2]*20 + [3]*20 + [4]*20)
    score = proteinnetworks.insight.getZScore(expected, generated)
    assert score > 0


def test_getzscore_bad_communities():

    # we expect most random guesses to place the community barrier closer to the we generated
    expected = np.asarray([1]*60 + [2]*40) # + [1]*20)
    generated = np.asarray([1]*5 + [2]*95) # + [3]*15 + [4]*35)
    score = proteinnetworks.insight.getZScore(expected, generated, numTrials=1000)
    assert score < 0


"""
Tests for getShannonEntropy

inputs -> a partition  (assuming labelled [1,m])
outputs -> the shannon entropy,
H(A) = - sum_i ( n_i/N * log(n_i/N).

Tests:
- invalid input (anything not an array of numbers ranging from 1 to m without gaps)
- all nodes in their own community
- all nodes in the same community
- "normal" community structure  
"""

def test_getshannonentropy_invalid_inputs():
    """Test that an array not containing [1,M] is rejected."""
    testarray = [-1]*5 + [2]*5 +[4]*5
    with pytest.raises(ValueError):
        H = proteinnetworks.insight.getShannonEntropy(testarray)

def test_getshannonentropy_single_community():
    """Test the limiting case of all nodes in one community."""
    testarray = [1]*100
    H = proteinnetworks.insight.getShannonEntropy(testarray)
    assert H==0

def test_getshannonentropy_unit_community():
    """Test the limiting case of all nodes in their own community."""
    testarray = range(1,101)
    H = proteinnetworks.insight.getShannonEntropy(testarray)
    assert math.isclose(H,-math.log(1/100))

def test_getshannonentropy_normal_community():
    """Test a traditional community structure."""
    testarray = [1]*10+[2]*10
    H = proteinnetworks.insight.getShannonEntropy(testarray)
    assert math.isclose(H,-math.log(1/2))



"""
tests for getMutualInfo

inputs -> two partitions (numpy arrays, labels in [1,m])
outputs -> mutual information, should always >=0 and < the average Shannon entropy.

tests:
- on inputs (check proper labelling)

- all nodes in same community
- all nodes in different communities
- matching communities
- completely disparate communities
"""
def test_getmutualinfo_invalid_inputs_labels():
    """Test that improperly labelled arrays are rejected.""" 
    testarray1 = [-1]*5 + [2]*5 +[4]*5
    testarray2 = [-1]*5 + [2]*5 +[4]*5

    with pytest.raises(ValueError):
        I = proteinnetworks.insight.getMutualInfo(testarray1, testarray2)

def test_getmutualinfo_invalid_inputs_lengths():
    """Test that different length arrays will be rejecteed."""
    testarray1 = [1]*5 + [2]*5 +[3]*5
    testarray2 = [1]*5 + [2]*5 

    with pytest.raises(TypeError):
        I = proteinnetworks.insight.getMutualInfo(testarray1, testarray2)

def test_getmutualinfo_single_community():
    """Test that single communities have I=0."""
    testarray1 = [1]*100
    testarray2 = [1]*100
    I = proteinnetworks.insight.getMutualInfo(testarray1, testarray2)
    assert I==0

def test_getmutualinfo_unit_community():
    """Test that unit communities have I= log(N)."""
    testarray1 = list(range(1,101))
    testarray2 = list(range(1,101))
    I = proteinnetworks.insight.getMutualInfo(testarray1, testarray2)
    assert math.isclose(I,math.log(100))


def test_getmutualinfo_matching_communities():
    """Test that matching bipartitions have I=log(N).""" 
    testarray1 = [1]*10+[2]*10
    testarray2 = [1]*10+[2]*10
    I = proteinnetworks.insight.getMutualInfo(testarray1, testarray2)
    assert math.isclose(I,math.log(2))


"""
tests for getNMI

input -> two partitions (invalid inputs will be caught by other functions)
output -> a number that should always be between 0 and 1.
"""

def test_getnmi_perfect_match():
    """Test that perfectly matching communities have NMI=1."""
    testarray1 = [1]*10+[2]*10
    testarray2 = [1]*10+[2]*10
    nmi = proteinnetworks.insight.getNMI(testarray1, testarray2)
    assert nmi == 1

def test_getnmi_perfect_match():
    """Test that random communities have NMI between 0 and 1."""
    numberOfTrials = 100
    for i in range(numberOfTrials):
        arraySize = np.random.randint(2,300)
        numberOfCommunites = np.random.randint(3,20)
        randomPartition1 = []
        randomPartition2 = []
        i=1
        while len(randomPartition1) < arraySize:
            communitySize = np.random.randint(1, arraySize/2+1)
            randomPartition1 += [i]*communitySize
            i+=1
        randomPartition1 = randomPartition1[:arraySize]
        i=1
        while len(randomPartition2) < arraySize:
            communitySize = np.random.randint(1, arraySize/2+1)
            randomPartition2 += [i]*communitySize
            i+=1
        randomPartition2 = randomPartition2[:arraySize]

        np.random.shuffle(randomPartition1)
        np.random.shuffle(randomPartition2)

        nmi = proteinnetworks.insight.getNMI(randomPartition1, randomPartition2)
        assert nmi >= 0 and nmi <= 1


"""
tests for getConductanceFromPartition

inputs -> one Network, one Partition
outputs -> a List of floats

tests:
invalid inputs:
    - partition and network are for different systems
    - other errors should be caught by subfunctions
- disconnected network
- fully connected network
- random partitions on the same network?
"""

def test_getconductancefrompartition_validinputs():
    """Test a success case where the network and partition correspond to the same system."""
    db = proteinnetworks.database.Database(password="bla")
    inputArgs = {
        "scaling": 4.5,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": "1ubq",
        "database": db
    }
    network = proteinnetworks.network.Network(**inputArgs)
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId(network.edgelistid),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    phis = proteinnetworks.insight.getConductanceFromPartition(network, partition)   
    assert phis

def test_getconductancefrompartition_invalidinputs():
    """Test a failure case where the network and partition correspond to different systems."""
    db = proteinnetworks.database.Database(password="bla")
    inputArgs = {
        "scaling": 4.5,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": "2vcr",
        "database": db
    }
    network = proteinnetworks.network.Network(**inputArgs)
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId("58dbe03fef677d54224a01da"),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    with pytest.raises(ValueError):
        phis = proteinnetworks.insight.getConductanceFromPartition(network, partition)   
    

"""
tests for getConductanceFromNodeSubset

inputs -> a 1d numpy array of true/false values and a 2d adjacency matrix

outputs -> a float in [0,1]

tests:
invalid inputs:
    - malformed node_subset?
    - malformed adjacency matrix?
- disconnected adj_matrix
- fully connected adj_matrix
- random node subsets 
"""

def test_getconductancefromnodesubset_incorrect_dimensions():
    """Check that the node subset only has indices that occur in the adjacency matrix."""
    adjacency  = np.asarray([[0,0,0,1],[0,0,1,1],[0,1,0,0],[1,1,0,0]])
    node_subset = [-1,0]
    with pytest.raises(ValueError):
        Q = proteinnetworks.insight.getConductanceFromNodeSubset(node_subset, adjacency)
    node_subset = [0,4]
    with pytest.raises(ValueError):
        Q = proteinnetworks.insight.getConductanceFromNodeSubset(node_subset, adjacency)

def test_getconductancefromnodesubset_non_numpy_array():
    adjacency  = [[0,0,0,1],[0,0,1,1],[0,1,0,0],[1,1,0,0]]
    node_subset = [0,1]
    with pytest.raises(TypeError):
        Q = proteinnetworks.insight.getConductanceFromNodeSubset(node_subset, adjacency)
    
def test_getconductancefromnodesubset_asymmetric_array():
    adjacency  = np.asarray([[0,0,0,1],[0,0,1,1],[0,1,0,0],[1,1,1,0]])
    node_subset = [0,1]
    with pytest.raises(ValueError):
        Q = proteinnetworks.insight.getConductanceFromNodeSubset(node_subset, adjacency)

def test_getconductancefromnodesubset_non_square_array():
    adjacency  = np.asarray([[0,0,0,1],[0,0,1,1],[0,1,0,0],[1,1,0,0], [0,0,0,0]])
    node_subset = [0,1]
    with pytest.raises(ValueError):
        Q = proteinnetworks.insight.getConductanceFromNodeSubset(node_subset, adjacency)

def test_getconductancefromnodesubset_negative_weights():
    adjacency  = np.asarray([[0,0,0,1],[0,0,1,1],[0,1,0,0],[1,1,0,-1]])
    node_subset = [0,1]
    with pytest.raises(ValueError):
        Q = proteinnetworks.insight.getConductanceFromNodeSubset(node_subset, adjacency)


def test_getconductancefromnodesubset_disconnected_graph():
    """Test that the conductance of a disconnected componenet is 0"""
    arrays = [np.ones((5,5)), np.ones((5,5))]
    adjacency = block_diag(*arrays)
    np.fill_diagonal(adjacency,0)
    node_subset = list(range(5))
    phi = proteinnetworks.insight.getConductanceFromNodeSubset(node_subset, adjacency)
    assert phi == 0

def test_getconductancefromnodesubset_fullyconnected_graph():
    """Test that the conductance of a fully-connected subgraph is 1"""
    arrays = [np.ones((5,5)), np.ones((5,5))]
    adjacency = block_diag(*arrays)
    # np.fill_diagonal(adjacency,0)
    adjacency = np.ones((10,10)) - adjacency
    node_subset = list(range(5))
    phi = proteinnetworks.insight.getConductanceFromNodeSubset(node_subset, adjacency)
    assert phi == 1

def test_getconductancefromnodesubset_random_graphs():
    """Test that for many random graphs and subgraphs the conductance lies within limits."""
    numTrials = 10
    for i in range(numTrials):
        l = np.random.randint(1,10)
        k = np.random.randint(10,50)
        p = np.random.random_sample()
        G = nx.relaxed_caveman_graph(l,k,p)
        # get the giant component
        adj = nx.convert_matrix.to_numpy_array(G)
        np.fill_diagonal(adj,0)
        node_subset = np.random.randint(0, G.number_of_nodes(), np.random.randint(1, G.number_of_nodes()))
        phi = proteinnetworks.insight.getConductanceFromNodeSubset(node_subset, adj)
        assert  phi >= 0




"""
tests for getModularityFromPartition.

for every level in the partition, generate the adjacency matrix then find the modularity.
Just test one successful path, others should be caught be subfunctions. 
"""

def test_getmodularityfrompartition_success():
    db = proteinnetworks.database.Database(password="bla")

    inputArgs = {
        "scaling": 4.5,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": "1ubq",
        "database": db
    }
    network = proteinnetworks.network.Network(**inputArgs)
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    Qs = proteinnetworks.insight.getModularityFromPartition(network,partition)
    assert len(Qs) == 2

"""
tests for getModularityFromAdjacencyMatrix.

Input -> an adjacency matrix (assumed to be a 2d numpy array)


Tests - a non-numpy array?
- An array with negative weights?
- A non-square array
- A disconnected graph
- A normal graph 
"""

def test_getmodularityfromadjacencymatrix_nonnumpy_array():
    """tests that anything other than a numpy array is rejected"""
    adjacency  = [[0,0,0,1],[0,0,1,1],[0,1,0,0],[1,1,0,0]]
    comlist = [1,1,0,0]
    with pytest.raises(TypeError):
        Q = proteinnetworks.insight.getModularityFromAdjacencyMatrix(adjacency, comlist)

def test_getmodularityfromadjacencymatrix_wrong_dimensions():
    """tests that if the comlist and adjacency matrix have different lengths an error is thrown"""
    adjacency  = np.asarray([[0,0,0,1],[0,0,1,1],[0,1,0,0],[1,1,0,0]])
    comlist = [1,1,0]
    with pytest.raises(ValueError):
        Q = proteinnetworks.insight.getModularityFromAdjacencyMatrix(adjacency, comlist)

def test_getmodularityfromadjacencymatrix_non_square_array():
    """tests that non-square arrays are rejected"""
    adjacency  = np.asarray([[0,0,0,1],[0,0,1,1],[0,1,0,0],[1,1,0,0], [1,1,0,0]])
    comlist = [1,1,0,0]
    with pytest.raises(ValueError):
        Q = proteinnetworks.insight.getModularityFromAdjacencyMatrix(adjacency, comlist)

def test_getmodularityfromadjacencymatrix_assymmetric_array():
    """tests that non-symmetric arrays are rejected"""
    adjacency  = np.asarray([[0,0,0,1],[0,0,1,1],[0,1,0,0],[1,1,1,0]])
    comlist = [1,1,0,0]
    with pytest.raises(ValueError):
        Q = proteinnetworks.insight.getModularityFromAdjacencyMatrix(adjacency, comlist)



def test_getmodularityfromadjacencymatrix_negative_weights():
    """tests that negative weights are rejected"""
    adjacency  = np.asarray([[0,0,0,1],[0,0,1,1],[0,1,0,0],[1,1,0,-1]])
    comlist = [1,1,0,0]
    with pytest.raises(ValueError):
        Q = proteinnetworks.insight.getModularityFromAdjacencyMatrix(adjacency, comlist)


def test_getmodularityfromadjacencymatrix_disconnected_graph():
    """test the modularity returned from two disconnected communities"""
    arrays = [np.ones((5,5)), np.ones((5,5))]
    adjacency = block_diag(*arrays)
    np.fill_diagonal(adjacency,0)
    comlist = [1]*5 + [2]*5
    Q = proteinnetworks.insight.getModularityFromAdjacencyMatrix(adjacency, comlist)
    assert math.isclose(Q,0.5)

    
def test_getmodularityfromadjacencymatrix_normal_communities():
    """test the modularity of random graphs, assert it is [-1, 1]"""
    numTrials = 100
    for i in range(numTrials):
        l = np.random.randint(1,10)
        k = np.random.randint(10,50)
        p = np.random.random_sample()
        G = nx.relaxed_caveman_graph(l,k,p)
        adj = nx.convert_matrix.to_numpy_array(G)
        comlist = []
        for i in range(l):
            comlist += [i]*k
        Q = proteinnetworks.insight.getModularityFromAdjacencyMatrix(adj, comlist)
        assert Q <= 1 and Q > -1


"""
tests for edgelistToGraph

input -> a weighted edgelist
output -> an undirected Graph

Tests:
- malformed edgelist (not rows of three)
- negative or non-float weights
-  i, j not integers
- normal edgelist
"""

def test_edgelisttograph_malformed_edgelist():
    """Tests that something other than rows of three throws an error."""
    edgelist = [[0,1,1], [1,2,4], [4,5]]
    with pytest.raises(ValueError):
        G = proteinnetworks.insight.edgelistToGraph(edgelist)

def test_edgelisttograph_incorrect_weights():
    """Tests that non-float or negative weights throw an error."""
    edgelist = [[0,1,"bla"], [1,2,4], [4,5,3]]
    with pytest.raises(ValueError):
        G = proteinnetworks.insight.edgelistToGraph(edgelist)
    edgelist = [[0,1,-1], [1,2,4], [4,5,3]]
    with pytest.raises(ValueError):
        G = proteinnetworks.insight.edgelistToGraph(edgelist)
   
def test_edgelisttograph_noninteger_nodelabels():
    """Tests that non-integer node labels will cause a ValueError."""
    edgelist = [["bla",1,1], [1,2,4], [4,5,3]]
    with pytest.raises(ValueError):
        G = proteinnetworks.insight.edgelistToGraph(edgelist)

def test_edgelisttograph_normal_edgelist():
    """Tests that a well-formatted edgelist correctly returns a Graph."""
    edgelist = [[1,2,1], [2,3,1], [3,1,1]]
    G = proteinnetworks.insight.edgelistToGraph(edgelist)
    assert G.number_of_edges() == 3
    assert G.number_of_nodes() == 3