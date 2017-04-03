"""Stores functionality related to the generation and analysis of community structures."""

import sys
import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt
from palettable.colorbrewer.qualitative import Set3_12
from .database import Database


class Partition:
    """
    Holds a community structure for a network (i.e the partition) and its parameters.

    Offers partition inspection and visualisation methods.
    """

    def __init__(self,
                 pdbref,
                 edgelistid,
                 detectionmethod,
                 r=-1,
                 N=-1,
                 database=None):
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

        if database:
            self.database = database
        else:
            try:
                self.database = Database()
                print("successfully connected")
            except IOError:
                print("Couldn't connect to server")
                sys.exit()
        # Attempt to extract the partition matching the given params
        doc = self.database.extractPartition(pdbref, edgelistid,
                                             detectionmethod, r, N)
        if doc:
            self.data = doc['data']
            self.partitionid = doc['_id']
            print("partition found")
        else:
            print("no partition fitting those parameters found: generating")

            data = self.generatePartition(pdbref, edgelistid, detectionmethod,
                                          r, N)
            self.data = data

            self.partitionid = self.database.depositPartition(
                pdbref, edgelistid, detectionmethod, r, N, data)

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
            str(N), "--silent"
        ])
        partition = treeFileToNestedLists("temp.tree")
        # Remove temporary files
        os.remove("temp.dat")
        os.remove("temp.tree")
        return partition

    def plotStripeDiagram(self, includePFAMDomains=False):
        """
        Plot the partition as a set of "stripe" plots.

        The includePFAMDomains flag plots the known PFAM structure alongside.
        """
        # Convert the list (or nested list) to a numpy array, for use with imshow.
        stripes = np.asarray(self.data, dtype=int)
        numPlots = np.shape(stripes)[0] + 1 if includePFAMDomains else np.shape(stripes)[0]

        fig, axes = plt.subplots(nrows=numPlots, figsize=(5, 5), sharex=True)

        if includePFAMDomains:
            pfamDomainArray = self.getPFAMDomainArray()
            assert len(pfamDomainArray) == np.shape(stripes)[1]
            stripes = np.vstack((pfamDomainArray, stripes))
            axes[0].imshow(
                np.vstack(
                    2 * [stripes[0, :]]),  # vstack otherwise imshow complains
                aspect=10,
                cmap='viridis')
            axes[0].xaxis.set_ticks_position('bottom')
            axes[0].yaxis.set_visible(False)
            axes[0].set_title("PFAM domain structure")
            axes[1].set_title("Generated structure")
            for i, ax in enumerate(axes[1:]):
                ax.imshow(
                    np.vstack(
                        2 * [stripes[i + 1, :]]),  # vstack otherwise imshow complains
                    aspect=10,
                    cmap=Set3_12.mpl_colormap,
                    interpolation="nearest")
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_visible(False)
        else:
            for i, ax in enumerate(axes):
                ax.imshow(
                    np.vstack(
                        2 * [stripes[i, :]]),  # vstack otherwise imshow complains
                    aspect=10,
                    cmap=Set3_12.mpl_colormap,
                    interpolation="nearest")
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_visible(False)

        plt.suptitle(self.pdbref)
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        plt.xlabel('Residue number')
        plt.show()

    def getPFAMDomainArray(self):
        """
        Get an array corresponding to the PFAM domains for a protein.

        Indexed by residue number, for now.
        """
        residues = []
        # Get the chain ID, start residue and end residue for the protein.
        mappings = self.database.extractMappings(self.pdbref, mappingtype="PFAM")
        if not mappings:
            raise ValueError("No PFAM data found for protein:", self.pdbref)
        else:
            residues = [x['data'] for x in mappings]

        # Read the pdb file, find the actual residue numbers.
        nodes = []
        pdb = self.database.extractPDBFile(self.pdbref)
        if not pdb:
            pdb = self.database.fetchPDBFileFromWeb(self.pdbref)
        for residue in residues:
            chain = residue['chainid']
            startResidue = int(residue['startresidue'])
            endResidue = int(residue['endresidue'])
            residueCounter = 0
            prevRes = 0
            firstNode = -1
            lastNode = -1
            for line in pdb:
                if line == "ENDMDL":
                    break
                if line[0:4] == "ATOM":
                    residueNumber = int(line[22:26])
                    if residueNumber != prevRes:
                        residueCounter += 1
                    prevRes = residueNumber
                    # now we have the node index.
                    if line[21] == chain and residueNumber == startResidue:
                        firstNode = residueCounter
                        startResidue = -1
                    if line[21] == chain and residueNumber == endResidue:
                        lastNode = residueCounter
                        break
            nodes.append([firstNode, lastNode])

        # Get the size of the array, given that the list may be nested
        n = len(self.data) if not any(isinstance(i, list) for i in self.data) else len(self.data[0])
        expectedDomains = np.ones(n, dtype=int)
        counter = 2
        for domain in nodes:
            startingindex = int(domain[0]) - 1
            endindex = int(domain[1])
            expectedDomains[startingindex:endindex] = counter
            counter += 1
        return expectedDomains


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
