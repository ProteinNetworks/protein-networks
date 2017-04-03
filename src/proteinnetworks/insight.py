"""
Functions for analysing the information content of the community structures.

Includes the SuperNetwork class, which is initialised by passing a Partition (which
contains a reference to the edgelist used).

The modified Jaccard for each level of the partition is calculated, and the level with the
best correspondence to the PFAM domains is chosen to generate a supernetwork.
"""
import numpy as np
import networkx as nx
import sys
import matplotlib.pyplot as plt
import warnings


class SuperNetwork:
    """
    A network generated from the community structure of the protein.

    Pull from the database if possible: otherwise generate anew.
    """

    def __init__(self, inputpartition):
        """Generate the network from an existing Partition."""
        # Get the input partition and edgelist
        self.pdbref = inputpartition.pdbref  # Save the details on the partition used
        self.database = inputpartition.database
        self.partitionid = inputpartition.partitionid

        # Attempt to extract the supernetwork matching the given params
        doc = self.database.extractSuperNetwork(self.pdbref, self.partitionid)

        if doc:
            self.data = doc['data']
            self.level = doc['level']
            print("supernetwork found")
        else:
            partition = inputpartition.data
            edgelist = inputpartition.database.extractDocumentGivenId(
                inputpartition.edgelistid)['data']

            # Find the level of the partition (assuming this is Infomap) with the best Jaccard
            try:
                pfamDomains = np.asarray(
                    inputpartition.getPFAMDomainArray(), dtype=int)
            except ValueError:
                print("No PFAM entry -> cannot generate supernetwork")
                raise ValueError

            maxJaccard = -1
            maxI = -1
            for i, col in enumerate(partition):
                jaccard = getModifiedJaccard(
                    pfamDomains, np.asarray(
                        col, dtype=int))
                print("Level {} has Jaccard {}".format(i, jaccard))
                if jaccard > maxJaccard:
                    maxJaccard = jaccard
                    maxI = i
            print("Using level {}".format(maxI))
            self.level = maxI
            partition = partition[maxI]

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
            self.database.depositSuperNetwork(self.pdbref, self.partitionid, self.level, self.data)

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
    numPFAMdomains = len(
        set(expectedArray))  # NB this include "1", the base counter
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
