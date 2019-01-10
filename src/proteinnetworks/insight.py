"""
Functions for analysing the information content of the community structures.

Includes the SuperNetwork class, which is initialised by passing a Partition (which
contains a reference to the edgelist used).

The modified Jaccard for each level of the partition is calculated, and the level with the
best correspondence to the PFAM domains is chosen to generate a supernetwork.
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import warnings
import itertools
import math
import logging
from .partition import Partition
from .network import Network
from typing import List

loggingLevels = {0: logging.ERROR, 1: logging.WARNING, 2: logging.INFO, 3: logging.DEBUG}


class SuperNetwork:
    """
    A network generated from the community structure of the protein.

    Pull from the database if possible: otherwise generate anew.
    Includes a single-chain reference, for use in the getIsomorphs function
    """

    def __init__(self, inputPartition, level=None, verbosity=1):
        """Generate the network from an existing Partition."""
        # Get the input partition and edgelist
        self.pdbref = inputPartition.pdbref  # Save the details on the partition used
        self.database = inputPartition.database
        self.partitionid = inputPartition.partitionid

        # Reset the verbosity
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        logger = logging.getLogger(__name__)
        logger.setLevel(loggingLevels[verbosity])
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

        ch = logging.StreamHandler()
        # ch.setFormatter()
        logger.addHandler(ch)

        partition = inputPartition.data
        edgelistDoc = inputPartition.database.extractDocumentGivenId(
            inputPartition.edgelistid)
        edgelist = edgelistDoc['data']
        self.chainref = edgelistDoc.get(
            'chainref')  # None if no chainref found
        # If no level is given, try to find the level best matching PFAM
        if level is None:
            try:
                pfamDomains = np.asarray(
                    inputPartition.getPFAMDomainArray(), dtype=int)
            except ValueError:
                logger.error("No PFAM entry -> cannot generate supernetwork")
                raise ValueError

            maxJaccard = -1
            maxI = -1
            for i, col in enumerate(partition):
                jaccard = getModifiedJaccard(pfamDomains,
                                             np.asarray(col, dtype=int))
                logger.info("Level {} has Jaccard {}".format(i, jaccard))
                if jaccard > maxJaccard:
                    maxJaccard = jaccard
                    maxI = i
            logger.info("Using level {}".format(maxI))
            self.level = maxI

        else:
            logger.info("Using specified level:", level)
            self.level = int(level)

        partition = partition[self.level]

        # Attempt to extract the supernetwork matching the given params
        doc = self.database.extractSuperNetwork(self.pdbref, self.partitionid,
                                                level)

        if doc:
            self.data = doc['data']
            logger.info("supernetwork found")

        else:
            # Generate the supernetwork
            communityEdgeList = {}
            for i, j, _ in edgelist:
                com_i, com_j = partition[int(i) - 1], partition[int(j) - 1]
                if com_i != com_j:
                    if not (com_i, com_j) in communityEdgeList:
                        if (com_j, com_i) in communityEdgeList:
                            communityEdgeList[(com_j, com_i)] += 1
                        else:
                            communityEdgeList[(com_i, com_j)] = 1
                    else:
                        communityEdgeList[(com_i, com_j)] += 1

            communityEdgeListSorted = []
            for row, weight in communityEdgeList.items():
                i, j = row
                communityEdgeListSorted.append([i, j, weight])
            communityEdgeListSorted.sort()

            self.data = communityEdgeListSorted
            self.database.depositSuperNetwork(self.pdbref, self.partitionid,
                                              self.level, self.data)

    @classmethod
    def fromPartitionId(SuperNetwork, partitionid, database, level=None):
        """
        Given a database and a partitionid, generate the Partition class.

        Then generate the SuperNetwork as normal from the partition.
        FIXME: this is really very convoluted.
        """
        partitionDetails = database.extractDocumentGivenId(partitionid)

        if 'r' in partitionDetails:
            inputPartition = Partition(
                partitionDetails['pdbref'],
                partitionDetails['edgelistid'],
                partitionDetails['detectionmethod'],
                r=partitionDetails['r'],
                database=database)
        elif 'N' in partitionDetails:
            inputPartition = Partition(
                partitionDetails['pdbref'],
                partitionDetails['edgelistid'],
                partitionDetails['detectionmethod'],
                N=partitionDetails['N'],
                database=database)
        return SuperNetwork(inputPartition=inputPartition, level=level)

    def draw(self):
        """Draw the reduced edgelist using NetworkX."""
        G = nx.Graph()
        for i, j, weight in self.data:
            G.add_edge(i, j, weight=weight)

        pos = nx.spring_layout(G, k=10)
        fig, ax = plt.subplots(figsize=(5, 5))

        # Suppress MPL's complaining, as it's a NetworkX problem.
        warnings.filterwarnings("ignore")
        nx.draw(G, pos=pos, node_color="grey")
        ax.set_title("Community network for {}".format(self.pdbref))
        plt.show()

    def getIsomorphs(self, subset=None):
        """Get all proteins in the database with an isomorphic supernetwork."""
        # Generate the NetworkX graph for the supernetwork
        G = nx.Graph()
        for i, j, weight in self.data:
            G.add_edge(i, j, weight=weight)

        # Get a cursor for all supernetworks in the database
        if subset is None:
            proteins = self.database.extractAllSuperNetworks(
                pdbref=self.pdbref)
        else:
            proteins = subset
        isomorphs = []
        for protein in proteins:
            G2 = nx.Graph()
            if type(protein) is dict:
                for i, j, weight in protein['data']:
                    G2.add_edge(i, j, weight=weight)
                if nx.faster_could_be_isomorphic(G, G2) and nx.is_isomorphic(
                        G, G2):
                    isomorphs.append(protein['pdbref'])
            elif type(protein) is SuperNetwork:
                for i, j, weight in protein.data:
                    G2.add_edge(i, j, weight=weight)
                if nx.faster_could_be_isomorphic(G, G2) and nx.is_isomorphic(
                        G, G2):
                    isomorphicProtein = protein.pdbref
                    if protein.chainref is not None:
                        isomorphicProtein += "_{}".format(protein.chainref)
                    isomorphs.append(isomorphicProtein)
        return isomorphs

    def getWeakIsomorphs(self, subset=None):
        """
        Get all proteins in the database with an weakly isomorphic supernetwork.

        Returns a list [self.pdbref, otherpdbref, simscore] for all proteins with
        a simscore > 0.5.

        If a subset of the supernetworks are given, this is used. (subset must be a numpy array)
        """
        # Generate the NetworkX graph for the supernetwork
        G = nx.Graph()
        for i, j, weight in self.data:
            G.add_edge(i, j)

        G = nx.convert_node_labels_to_integers(G)
        # Get a cursor for all supernetworks in the database

        if subset.any():
            proteins = subset
        else:
            proteins = self.database.extractAllSuperNetworks(
                pdbref=self.pdbref)
        weakIsomorphs = []
        for protein in proteins:
            G2 = nx.Graph()
            for i, j, weight in protein['data']:
                G2.add_edge(i, j)

            G2 = nx.convert_node_labels_to_integers(G2)
            # Get the maximum common subgraph for the two supernetworks
            try:
                MCS = getMCS(G, G2)
            except ValueError:
                continue
            if MCS:
                similarity = MCS.number_of_nodes() / (
                    max(G.number_of_nodes(), G2.number_of_nodes()))
            else:
                similarity = 0

            if similarity > 0.5:
                weakIsomorphs.append(
                    [self.pdbref, protein['pdbref'],
                     str(similarity)])

        return weakIsomorphs


class SuperNetworkNullModel:
    """
    A network generated from a null model, using a given community structure as a base.

    Never store in the database
    Includes a single-chain reference, for use in the getIsomorphs function
    """

    def __init__(self, inputPartition, level=None):
        """Generate the network from an existing Partition."""
        # Get the input partition and edgelist
        self.pdbref = inputPartition.pdbref  # Save the details on the partition used
        self.database = inputPartition.database
        self.partitionid = inputPartition.partitionid
        partition = inputPartition.data
        edgelistDoc = inputPartition.database.extractDocumentGivenId(
            inputPartition.edgelistid)
        edgelist = edgelistDoc['data']
        self.chainref = edgelistDoc.get(
            'chainref')  # None if no chainref found
        # If no level is given, try to find the level best matching PFAM

        partition = generateNullModel(partition[0])
        # Generate the supernetwork
        communityEdgeList = {}
        for i, j, _ in edgelist:
            com_i, com_j = partition[int(i) - 1], partition[int(j) - 1]
            if com_i != com_j:
                if not (com_i, com_j) in communityEdgeList:
                    if (com_j, com_i) in communityEdgeList:
                        communityEdgeList[(com_j, com_i)] += 1
                    else:
                        communityEdgeList[(com_i, com_j)] = 1
                else:
                    communityEdgeList[(com_i, com_j)] += 1

        communityEdgeListSorted = []
        for row, weight in communityEdgeList.items():
            i, j = row
            communityEdgeListSorted.append([i, j, weight])
        communityEdgeListSorted.sort()

        self.data = communityEdgeListSorted

    def draw(self):
        """Draw the reduced edgelist using NetworkX."""
        G = nx.Graph()
        for i, j, weight in self.data:
            G.add_edge(i, j, weight=weight)

        pos = nx.spring_layout(G, k=10)
        fig, ax = plt.subplots(figsize=(5, 5))

        # Suppress MPL's complaining, as it's a NetworkX problem.
        warnings.filterwarnings("ignore")
        nx.draw(G, pos=pos, node_color="grey")
        ax.set_title("Community network for {}".format(self.pdbref))
        plt.show()

    def getIsomorphs(self, subset=None):
        """Get all proteins in the database with an isomorphic supernetwork."""
        # Generate the NetworkX graph for the supernetwork
        G = nx.Graph()
        for i, j, weight in self.data:
            G.add_edge(i, j, weight=weight)
        # Get a cursor for all supernetworks in the database
        if subset is None:
            proteins = self.database.extractAllSuperNetworks(
                pdbref=self.pdbref)
        else:
            proteins = subset
        isomorphs = []
        for protein in proteins:

            G2 = nx.Graph()
            if type(protein) is dict:
                for i, j, weight in protein['data']:
                    G2.add_edge(i, j, weight=weight)
                if nx.faster_could_be_isomorphic(G, G2) and nx.is_isomorphic(
                        G, G2):
                    isomorphs.append(protein['pdbref'])
            elif type(protein) is SuperNetworkNullModel:
                for i, j, weight in protein.data:
                    G2.add_edge(i, j, weight=weight)
                if nx.faster_could_be_isomorphic(G, G2) and nx.is_isomorphic(
                        G, G2):
                    isomorphicProtein = protein.pdbref
                    if protein.chainref is not None:
                        isomorphicProtein += "_{}".format(protein.chainref)
                    isomorphs.append(isomorphicProtein)
        return isomorphs

    def getWeakIsomorphs(self, subset=None):
        """
        Get all proteins in the database with an weakly isomorphic supernetwork.

        Returns a list [self.pdbref, otherpdbref, simscore] for all proteins with
        a simscore > 0.5.

        If a subset of the supernetworks are given, this is used. (subset must be a numpy array)
        """
        # Generate the NetworkX graph for the supernetwork
        G = nx.Graph()
        for i, j, weight in self.data:
            G.add_edge(i, j)

        G = nx.convert_node_labels_to_integers(G)
        # Get a cursor for all supernetworks in the database

        if subset.any():
            proteins = subset
        else:
            proteins = self.database.extractAllSuperNetworks(
                pdbref=self.pdbref)
        weakIsomorphs = []
        for protein in proteins:
            G2 = nx.Graph()
            for i, j, weight in protein['data']:
                G2.add_edge(i, j)

            G2 = nx.convert_node_labels_to_integers(G2)
            # Get the maximum common subgraph for the two supernetworks
            try:
                MCS = getMCS(G, G2)
            except ValueError:
                continue
            if MCS:
                similarity = MCS.number_of_nodes() / (
                    max(G.number_of_nodes(), G2.number_of_nodes()))
            else:
                similarity = 0

            if similarity > 0.5:
                weakIsomorphs.append(
                    [self.pdbref, protein['pdbref'],
                     str(similarity)])

        return weakIsomorphs


def getModifiedJaccard(expectedArray, generatedArray):
    """
    A scoring function for each PFAM domain in a protein.

    Requires:
    - The PDB file
    - The PFAM/PDB mapping
    - The .tree file (or other community structure)

    Scoring algorithm:

    For each domain:
        - Collect all the generated modules that overlap with the PFAM domain
        - Calculate the Jaccard index:
            J = | A ∩ B | / | A ∪ B |
        - Return the mean Jaccard for all generated modules, weighted by the intersection.

    The final score is the mean value over all domains.
    """
    if not type(expectedArray) == np.ndarray or not type(generatedArray) == np.ndarray:
        raise TypeError("both inputs must be numpy arrays") 
    if not np.equal(np.mod(expectedArray, 1),0).all() or not np.equal(np.mod(generatedArray, 1),0).all():
        raise TypeError("both inputs must be arrays of integers") 
    if len(expectedArray) == 0 or len(generatedArray) == 0:
        raise ValueError("both arrays must have non-zero length")
    if len(expectedArray) != len(generatedArray):
        raise ValueError("both arrays must have the same length")

    numPFAMdomains = len(
        set(expectedArray))  # NB this include "1", the base counter
    if numPFAMdomains == 1:
        raise ValueError("Input expected array has no domains")
    jaccards = []
    for i in range(2, numPFAMdomains + 1):
        # Get the modules with some overlap.
        jaccard = []
        intersections = []
        overlappingModules = set(generatedArray[expectedArray == i])
        for j in overlappingModules:
            intersection = np.count_nonzero(
                np.logical_and(generatedArray == j, expectedArray == i))
            union = np.count_nonzero(
                np.logical_or(generatedArray == j, expectedArray == i))
            jaccard.append(intersection / union)
            intersections.append(intersection)

        # weight the terms according to the overlap proportion.
        jaccard = [
            x * y / sum(intersections) for x, y in zip(jaccard, intersections)
        ]
        jaccard = sum(jaccard)
        jaccards.append(jaccard)
    return np.mean(jaccards)


def getZScore(expectedArray, generatedArray, numTrials=100):
    """
    Get the z-score.

    Defined as: z = (J - mu ) / sigma . Find mu and sigma by generating
    null models.
    """

    if numTrials < 1:
        raise ValueError("number of trials must be a positive integer")

    J = getModifiedJaccard(expectedArray, generatedArray)

    nullJaccard = []
    for i in range(numTrials):
        nullmodel = generateNullModel(generatedArray)
        nullJaccard.append(getModifiedJaccard(expectedArray, nullmodel))

    mu = np.mean(nullJaccard)
    sigma = np.std(nullJaccard)
    if math.isclose(mu, J):
        return 0.0
    else:
        return (J - mu) / sigma


def generateNullModel(testPartition):
    """
    From a given partition, generate a null model.

    Here the null model has the same number of boundaries as the generated partition, but with
    the boundaries arbitrarily placed (and the same number of communities in total.)
    """
    # Get the total number of communities in the partition

    if not type(testPartition) == np.ndarray:
        raise TypeError("input must be a numpy array") 
    if not np.equal(np.mod(testPartition, 1),0).all():
        raise TypeError("input must be an array of integers") 
    if len(testPartition) == 0:
        raise ValueError("array must have non-zero length")

    numCommunities = len(set(testPartition))

    if set(testPartition) != set(range(1,numCommunities+1)):
        raise ValueError("communities must be labelled [1,m] where m is the number of communities")

    # Get the number of boundaries
    prevI = -1
    numBoundaries = -1  # Start from -1 to avoid counting the start as a boundary
    for i in testPartition:
        if i != prevI:
            numBoundaries += 1
        prevI = i

    # Place the boundaries arbitrarily, and assign each segment to a randomly chosen community
    newBoundaries = []
    for i in range(numBoundaries):
        randint = np.random.randint(low=1, high=len(testPartition))
        while randint in newBoundaries:
            randint = np.random.randint(low=1, high=len(testPartition))
        newBoundaries.append(randint)
    newBoundaries.sort()
    newBoundaries = [0] + newBoundaries + [len(testPartition) - 1]
    # Fill the gaps with communities in the range 1... numComs.

    # We want to randomly assign, but ensure the same number of communities.
    # So assign minimum number of communities, then add to reach quota.
    # FIXME this currently sometimes hangs forever, if e.g [1, 2, 1, 1]
    # FIXME then no way of subsequently shuffling to make it work

    newCommunities = [i + 1 for i in range(numCommunities)]
    while len(newCommunities) != len(newBoundaries) - 1:
        newCommunities.append(np.random.randint(numCommunities) + 1)

    # Shuffle until no two numbers are together
    np.random.shuffle(newCommunities)
    valid = False
    counter = 0
    while not valid:
        counter += 1
        if counter == 10:
            # FIXME dirty hack to break out of "too-many-repeats" problem case.
            newCommunities = [i + 1 for i in range(numCommunities)]
            while len(newCommunities) != len(newBoundaries) - 1:
                newCommunities.append(np.random.randint(numCommunities) + 1)
            counter = 0
        valid = True
        np.random.shuffle(newCommunities)
        for i in range(len(newCommunities) - 1):
            if newCommunities[i] == newCommunities[i + 1]:
                valid = False

    nullModel = np.zeros(len(testPartition), dtype=int)
    for i in range(len(newBoundaries) - 1):
        nullModel[newBoundaries[i]:newBoundaries[i + 1] + 1] = newCommunities[
            i]

    # Check the null model
    # Get the number of boundaries
    prevI = -1
    nullModelNumBoundaries = -1  # Start from -1 to avoid counting the start as a boundary
    for i in nullModel:
        if i != prevI:
            nullModelNumBoundaries += 1
        prevI = i
    assert set(nullModel) == set(testPartition)
    assert len(nullModel) == len(testPartition)
    if nullModelNumBoundaries != numBoundaries:
        logging.error(newBoundaries)
        logging.error(newCommunities)
        logging.error(testPartition.tolist())
        logging.error()
        logging.error(nullModel.tolist())
        raise NotImplementedError
    return nullModel


def getMCS(G1, G2):
    """Take two networkx graphs, return the MCS as a networkx graph."""

    if not type(G1) == nx.Graph or not type(G2) == nx.Graph:
        raise TypeError("inputs for graph isomorphism must be undirected graphs")

    # Let G1 be the smaller graph
    if G1.number_of_nodes() > G2.number_of_nodes():
        temp = G2
        G2 = G1
        G1 = temp

    G1 = nx.relabel.convert_node_labels_to_integers(G1)
    G2 = nx.relabel.convert_node_labels_to_integers(G2)

    N_G1 = G1.number_of_nodes()
    N_G2 = G2.number_of_nodes()
    if N_G2 > 35:
        raise ValueError("Graph is too large")
    nodelist_G1 = list(range(N_G1))
    nodelist_G2 = list(range(N_G2))
    for i in range(N_G1, 0, -1):
        # Get all choose(N_G1, i) possible selections of [1... N_G1]
        # print(i)
        for subgraph_G1_nodelist in itertools.combinations(nodelist_G1, i):
            subgraph_G1 = G1.subgraph(subgraph_G1_nodelist)
            # Check whether subgraph_G1 is isomorphic to any subgraph of the same size in G2
            for subgraph_G2_nodelist in itertools.combinations(nodelist_G2, i):
                subgraph_G2 = G2.subgraph(subgraph_G2_nodelist)
                if nx.is_isomorphic(subgraph_G1, subgraph_G2):
                    return nx.Graph(subgraph_G1)


def getShannonEntropy(partition):
    """
    H(A) = - sum_i ( n_i/N * log(n_i/N).

    Here n_i is the number of nodes in community i, and N the total number of nodes.
    """
    N = len(partition)
    H = 0
    totalComNumber = len(set(partition))
    comSize = np.zeros(totalComNumber, dtype=int)
    for i in range(1, totalComNumber + 1):
        for j in range(0, N):
            if partition[j] == i:
                comSize[i - 1] += 1
    logComSize = np.log(comSize / N)
    for i in range(0, totalComNumber):
        H -= comSize[i] / N * logComSize[i]
    return H


def getMutualInfo(generated, expected):
    """
    Get the Mutual Information.

    I(A,B) = = sum_i sum_j n_ij / N * log ( n_ij N / n_i n_j ).
    Here n_ij is the number of nodes of community i in partition A that are
    in community j of partition B.

    Refs:   Ronhovde, P. & Nussinov, Z.
                Multiresolution community detection for megascale networks
                by information-based replica correlations.
                Phys. Rev. E 80, 016109 (2009)
    """
    N1 = len(generated)
    N2 = len(expected)
    assert N1 == N2
    comNumber1 = len(set(generated))
    comNumber2 = len(set(expected))
    confusionmatrix = np.zeros((comNumber1, comNumber2), dtype=int)
    for i in range(N1):
        confusionmatrix[generated[i] - 1, expected[i] - 1] += 1
    rowsums = np.sum(confusionmatrix, axis=1)
    colsums = np.sum(confusionmatrix, axis=0)

    mutualinfo = 0
    for i in range(0, comNumber1):
        for j in range(0, comNumber2):
            if confusionmatrix[i, j] != 0:
                mutualinfo += confusionmatrix[i, j] / N1 * math.log(
                    confusionmatrix[i, j] * N1 / (rowsums[i] * colsums[j]))
    return mutualinfo


def getNMI(partitionA, partitionB):
    """Given two partitions, return the NMI."""
    HA = getShannonEntropy(partitionA)
    mutualinfo = getMutualInfo(partitionA, partitionB)
    HB = getShannonEntropy(partitionB)
    NMI = 2 * mutualinfo / (HA + HB)
    return NMI


def getConductanceFromPartition(network: Network,
                                partition: Partition) -> List[float]:
    """Given a Network and Partition, return a conductance for each level of the partition."""
    generatedArray = partition.data
    adjacency_matrix = network.getAdjacencyMatrix()

    # NB this assumes labelling from 1 to m
    conductances = []
    for col in generatedArray:
        conductance = [
            getConductanceFromNodeSubset(
                np.where(np.asarray(col, dtype=int) == j + 1)[0],
                adjacency_matrix) for j in range(len(set(col)))
        ]
        conductances.append(conductance)
    return conductances


def getConductanceFromNodeSubset(node_subset, adjacency_matrix):
    """
    Return conductance given an adjacency matrix and a node_subset.

    The node-subset is a 1D array with same dimension as the
    adjacency matrix.
    Conductance is defined for a subset S, and its complement Sbar :

    C = sum_{i in S, j in Sbar} (a_{ij}) / min (a(S), a(Sbar))

    where a is the adjacency matrix, and a(S) is sum_{i in S, j in V} a_ij
    i.e. the total weight of edges indicent with S.
    """
    node_complement = [
        i for i in range(len(adjacency_matrix)) if i not in node_subset
    ]
    numerator = 0
    for i in node_subset:
        for j in node_complement:
            numerator += adjacency_matrix[i, j]

    denominator_1 = 0
    for i in node_subset:
        for j in range(len(adjacency_matrix)):
            denominator_1 += adjacency_matrix[i, j]

    denominator_2 = 0
    for i in node_complement:
        for j in range(len(adjacency_matrix)):
            denominator_2 += adjacency_matrix[i, j]

    conductance = numerator / min(denominator_1, denominator_2)

    return conductance


def getModularityFromPartition(network: Network,
                               partition: Partition) -> List[float]:
    """Given a Network and Partition, return Newman's modularity for each level of the partition."""
    generatedArray = np.atleast_2d(np.asarray(partition.data, dtype=int))
    adjacency_matrix = network.getAdjacencyMatrix()
    # Calculate Q for each level
    Qs = []
    for col in generatedArray:

        Q = getModularityFromAdjacencyMatrix(adjacency_matrix, col)
        Qs.append(Q)
    return Qs


def getModularityFromAdjacencyMatrix(adj, comlist):
    r"""
    Given an adjacency matrix and a list of communities, return the modularity.

    Q = 1/2w \sum_{i,j}( A_ij - A_i A _j / 2w  \delta( C_i, C_j))

      = \sum_{s} (w_{ss}/w = (w_s/2w)**2)
    """
    rowsums = np.sum(adj, axis=0)
    N = len(adj)
    inversestrength = 1 / np.sum(adj)
    Q = 0.0
    for i in range(0, N):
        wi = rowsums[i]
        pi = comlist[i]
        for j in range(0, i):
            wj = rowsums[j]
            if pi == comlist[j]:  # if nodes i and j are in the same community
                Q += (adj[i, j] - wi * wj * inversestrength)

    ondiagonals = sum(
        [adj[i, i] - (rowsums[i]**2) * inversestrength for i in range(N)])

    return (Q * 2 + ondiagonals) * inversestrength


def edgelistToGraph(edgelist):
    """
    Take a weighted edgelist of the form [[sourcenode, targetnode, weight], ...]
    and return a networkx Graph.
    """
    G = nx.Graph()
    for i, j, weight in edgelist:
        G.add_edge(i, j, weight=weight)
    return G


def drawTree(tree):
    """
    Draw a Graph object using graphviz
    """
    # from networkx.drawing.nx_agraph import graphviz_layout
    warnings.filterwarnings("ignore")  # Until Networkx gets its act together

    # pos = graphviz_layout(tree, prog='twopi', root="0")
    plt.figure(figsize=(8, 8))
    nx.draw(
        tree,
        # pos,
        node_size=20,
        alpha=0.5,
        node_color="blue",
        with_labels=True)
    plt.show()


def getTreeSimilarity(partition1, partition2):
    """
    Given two Partitions, convert them into trees, then perform iterative labelling.

    This should reveal the similarity of the two trees.

    Algorithm:

    Starting at the penultimate level, label the subtrees by isomorphism class.
    Then do a labelled isomorphism check of the higher nodes, all the way to the top.
    """
    import json  # for pretty printing
    tree1 = toTree(partition1.data)
    # tree1 = edgelistToGraph(treeEdgelist1)
    # tree2 = toTree(partition2.data)
    # tree2 = edgelistToGraph(treeEdgelist2)

    # Get the leaves
    leafNodes1 = [x for x, y in tree1.degree() if y == 1]
    # leafNodes2 = [x for x, y in tree2.degree() if y == 1]
    # print(leafNodes1)
    print(json.dumps(list(tree1.nodes(data=True)), indent=2))

    # sizeOfTree = tree1.number_of_nodes()
    # numberOfLabelledNodes = 0
    for node in leafNodes1:
        tree1.nodes[node]['childrenLabels'] = []
        tree1.nodes[node]['isoLabel'] = 0
        tree1.nodes[node][
            "size"] = 1  # Assumes all initial leaf nodes are of size 1
        # numberOfLabelledNodes += 1
        bottomLevel = tree1.nodes[node]["level"]
    print(json.dumps(list(tree1.nodes(data=True)), indent=2))
    # print(sizeOfTree)
    # print(numberOfLabelledNodes)
    print(bottomLevel)
    for level in range(bottomLevel - 1, 0, -1):
        nodesAtLevel = [
            x for x, y in tree1.nodes(data=True) if y['level'] == level
        ]
        print(nodesAtLevel)
        for node in nodesAtLevel:
            print(tree1.nodes(data=True)[node])
            print(tree1.nodes[node]['children'])
    # while numberOfLabelledNodes < sizeOfTree:


def toTree(data):
    """
    5/2 OUTPUTS A NETWORKX GRAPH
    Take the partition, and output an edgelist for a  tree in which each node is a community,
    with edge weight the community size, connected to its parent node.
    For each community, index the node, the level, and the region of the network it spans.
    Then for each level, if a community in the previous level has an overlap, then create a
    link.
    Note we expect this to produce a simple tree i.e. one and only one community in the
    previous level should show this. SHOULD TEST
    """
    nodes = [{"label": 0, "level": 0, "region": [True] * len(data[0])}]

    G = nx.Graph()

    G.add_node(0, level=0, children=[])  # , region=[True] * len(data[0]))

    nodeLabelCounter = 1

    treeDepth = len(data)
    for i in range(treeDepth):
        partition = np.asarray(data[i], dtype=int)
        for community in set(partition):
            node = {
                "label": nodeLabelCounter,
                "level": i + 1,
                "region": list(partition == community)
            }
            G.add_node(nodeLabelCounter, level=i + 1, children=[])
            # , region=list(partition == community)) leave off for now
            nodes.append(node)
            nodeLabelCounter += 1

    # Now we have a list of "node" dicts. For each node, check whether any nodes in the previous
    # level overlap. If so, make an edge of weight = size of overlap.
    # NB this is a super inefficient way to do this
    edges = []
    for node in nodes:
        prevLevel = [x for x in nodes if x['level'] == node['level'] - 1]
        if not prevLevel:
            # Will occur for the root node
            continue
        for prevLevelNode in prevLevel:
            overlap = [
                x and y
                for x, y in zip(node['region'], prevLevelNode['region'])
            ]
            overlapSize = sum(overlap)
            if overlapSize != 0:
                edges.append(
                    [node['label'], prevLevelNode['label'], overlapSize])

                G.add_edge(
                    node['label'], prevLevelNode['label'], weight=overlapSize)
                # print(G.nodes[prevLevelNode['label']])
                # print(G.nodes[prevLevelNode['label']]['children'])
                # print(node['label'])
                G.nodes[prevLevelNode['label']]['children'].append(
                    node['label'])

    # return edges
    return G
