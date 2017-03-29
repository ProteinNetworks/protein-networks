"""Stores functionality related to the generation and analysis of community structures."""

import sys
import subprocess
import numpy as np
import os
from .database import Database


class Partition:
    """
    Holds a community structure for a network (i.e the partition) and its parameters.

    Offers partition inspection and visualisation methods.
    """

    def __init__(self, pdbref, edgelistid, detectionmethod, r=-1, N=-1):
        """
        Initialise the Partition with a given set of params.

        Pull from the database if possible, else generate anew. (should perhaps rethink
        this, since community detection is a slow process)
        Partition is specified by:
            - doctype == partition
            - edgelistid: The _id of the edgelist used in generating the partition
               ^ (should check this is valid)
            - detectionmethod: (AFG | Infomap) as it stands
            if detectionmethod == AFG:
                r : The AFG parameter used in generating the partition (I should refactor this)
            - PDB reference
        """
        self.pdbref = pdbref
        self.edgelistid = edgelistid
        self.detectionmethod = detectionmethod
        if r != -1:
            self.r = r
        if N != -1:
            self.N = N
        # Try to connect to the database
        try:
            self.database = Database()
            print("successfully connected")
        except IOError:
            print("Couldn't connect to server")
            sys.exit()
        # Attempt to extract the edgelist matching the given params
        doc = self.database.extractPartition(pdbref, edgelistid,
                                             detectionmethod, r, N)
        if doc:
            self.partition = doc['data']
            print("partition found")
        else:
            print("no partition fitting those parameters found: generating")

            data = self.generatePartition(pdbref, edgelistid, detectionmethod,
                                          r, N)
            self.data = data
            self.database.depositPartition(pdbref, edgelistid, detectionmethod,
                                           r, N, data)

    def generatePartition(self, pdbref, edgelistid, detectionmethod, r, N):
        """Generate a community structure using the parameters supplied."""
        # Get the network.
        edgelist = self.database.extractDocumentGivenId(edgelistid)['data']
        assert detectionmethod == "Infomap"  # for now
        assert N > 0
        # Write the edgelist to a temporary file.
        with open("temp.dat", mode='w') as flines:
            flines.write("\n".join(" ".join(map(str, x)) for x in edgelist))
        # Run Infomap on the edgelist
        subprocess.run([
            "Infomap", "temp.dat", ".", "-i", "link-list", "--tree", "-N",
            str(N)
        ])
        partition = treeFileToNestedLists("temp.tree")
        # Remove temporary files
        os.remove("temp.dat")
        os.remove("temp.tree")
        return partition


def treeFileToNestedLists(inputTreeFile):
    """
    Take a path to a .tree file, output a list of partitions.

    Input: A path to a .tree file (can be jagged)
    1:1:1:1
    1:1:1:2
    1:1:2:1
    1:1:2:2
    1:1:2:3
    ...

    Output:
    A np array (as list of lists) sorted according to the node index (last column is the node index)
    """
    inputArray = []
    maxLevels = 0
    with open(inputTreeFile, mode='r') as readfile:
        for line in readfile:
            if line[0] == '#':
                continue
            cols = line.split(" ")
            nodeindex = int(cols[-1])
            trees = [int(x) for x in cols[0].split(":")]
            numLevels = len(trees)
            if numLevels > maxLevels:
                maxLevels = numLevels
            trees.append(nodeindex)
            inputArray.append(trees)

    npArray = np.ones((len(inputArray), maxLevels + 1), dtype=int)

    for i, row in enumerate(inputArray):
        for j, element in enumerate(row[:-1]):
            npArray[i, j] = element
        npArray[i, -1] = row[-1]

    inputArray = npArray
    numNodes = inputArray.shape[0]
    numLevels = inputArray.shape[1] - 1  # Last column is the node index
    # Relabel the communities
    for i in range(1, numLevels):  # iterate over all sublevels
        prevCurrentLevelElement = 1
        prevSuperLevelElement = 1
        offset = 0
        for j in range(numNodes):
            currentLevelElement = inputArray[j, i]
            superLevelElement = inputArray[j, i - 1]
            if superLevelElement != prevSuperLevelElement:
                offset += prevCurrentLevelElement
            inputArray[j, i] += offset
            prevCurrentLevelElement = currentLevelElement
            prevSuperLevelElement = superLevelElement

    # Drop the last column, it's simply the node indices
    return inputArray[inputArray[:, -1].argsort()].T.tolist()[:-1]
