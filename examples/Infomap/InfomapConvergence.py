# coding: utf-8

# Check how the number of outer loops (N) changes the conductance
#  (as a proxy for "goodness") of the modules, for a variety of protein shapes and sizes
from bson.objectid import ObjectId
from typing import List
import numpy as np
import proteinnetworks


def getConductance(
        network: proteinnetworks.network.Network,
        partition: proteinnetworks.partition.Partition) -> List[float]:
    """Given a Network and Partition, return a conductance for each level of the partition."""
    generatedArray = partition.data
    adjacency_matrix = network.getAdjacencyMatrix()

    # NB this assumes labelling from 1 to m
    conductances = []
    for col in generatedArray:
        # Got to cast to a numpy array otherwise np.where falls over quietly
        conductance = [
            calculateConductance(np.where(np.asarray(col, dtype=int) == j + 1)[0], adjacency_matrix)
            for j in range(len(set(col)))
        ]
        conductances.append(conductance)
    return conductances


def calculateConductance(node_subset, adjacency_matrix):
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


db = proteinnetworks.database.Database(password="8moMlmkRnZurXsl")

# Test Proteins:
# 1fl3 1r76 2ghw 3mym 3rj3 3rpi

pdbRef = "1fl3"
inputArgs = {
    "scaling": 4.0,
    "edgelisttype": "residue",
    "hydrogenstatus": "noH",
    "pdbref": pdbRef,
    "database": db
}
proteinNetwork = proteinnetworks.network.Network(**inputArgs)

for N in [1, 10, 100, 1000, 10000]:
    partitionArgs = {
        "pdbref": pdbRef,
        "edgelistid": ObjectId(proteinNetwork.edgelistid),
        "detectionmethod": "Infomap",
        "N": N,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    getConductance(network=proteinNetwork, partition=partition)
