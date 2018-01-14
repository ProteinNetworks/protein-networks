"""Stores functionality related to the generation and analysis of community structures."""

import sys
import subprocess
import numpy as np
import os
import networkx as nx
import warnings
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
        numPlots = np.shape(
            stripes)[0] + 1 if includePFAMDomains else np.shape(stripes)[0]
        fig, axes = plt.subplots(nrows=numPlots, figsize=(10, 5), sharex=True)

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
                    np.vstack(2 * [stripes[i + 1, :]
                                   ]),  # vstack otherwise imshow complains
                    aspect=10,
                    cmap=Set3_12.mpl_colormap,
                    interpolation="nearest")
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_visible(False)
        else:
            for i, ax in enumerate(axes):
                ax.imshow(
                    np.vstack(2 * [stripes[i, :]
                                   ]),  # vstack otherwise imshow complains
                    aspect=10,
                    cmap=Set3_12.mpl_colormap,
                    interpolation="nearest")
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_visible(False)

        plt.suptitle(self.pdbref)
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        plt.xlabel('Residue number')
        plt.savefig("{}.pdf".format(self.pdbref), dpi=300)
        plt.show()

    def getPFAMDomainArray(self):
        """
        Get an array corresponding to the PFAM domains for a protein.

        Indexed by residue number, for now.
        """
        residues = []
        # Get the chain ID, start residue and end residue for the protein.
        mappings = self.database.extractMappings(
            self.pdbref, mappingtype="PFAM")
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
        n = len(self.data) if not any(isinstance(
            i, list) for i in self.data) else len(self.data[0])
        expectedDomains = np.ones(n, dtype=int)
        counter = 2
        for domain in nodes:
            startingindex = int(domain[0]) - 1
            endindex = int(domain[1])
            expectedDomains[startingindex:endindex] = counter
            counter += 1
        return expectedDomains

    def plotPymolStructure(self, level=-1, outputPng=False):
        """Plot the community structure overlaid onto the protein using PyMol.

        Generate the sequence of b-factor alterations for each level.
        Each column of the input Array alters a Pymol Object of the form
        $PDBREF_$INDEX where the pdbref is the 4-character identifier, and index the
        column number (starting from zeros)
        """
        pymolCommands = []

        # Work out whether the given edgelist is contact or atomic
        edgelisttype = self.database.extractDocumentGivenId(
            self.edgelistid)['edgelisttype']
        if edgelisttype == "residue":
            selector = "resi"
        elif edgelisttype == "atomic":
            selector = "index"

        # If a level if specified, plot only that level

        if level != -1:
            i = level
            col = self.data[i]
            numberOfCommunities = len(set(col))
            col = np.asarray(col, dtype=int)  # necessary for np.where()
            pymolCommand = "\n".join([
                "alter {0}_{1} and ({2} {3} ), b={4}".format(
                    self.pdbref, i, selector, " or {} ".format(selector).join(
                        str(x + 1) for x in np.where(col == com)[0]),
                    com / numberOfCommunities) for com in set(col)
            ])
            pymolCommands.append(pymolCommand)
        else:
            for i, col in enumerate(self.data):
                numberOfCommunities = len(set(col))
                col = np.asarray(col, dtype=int)  # necessary for np.where()
                pymolCommand = "\n".join([
                    "alter {0}_{1} and ({2} {3} ), b={4}".format(
                        self.pdbref, i, selector, " or {} ".format(selector).join(
                            str(x + 1) for x in np.where(col == com)[0]),
                        com / numberOfCommunities) for com in set(col)
                ])
                pymolCommands.append(pymolCommand)
        """
        Amalgamate the pymol commands so as to load a pdb file for each level of hierarchy,
        then run the bfactor alterations.
        """
        # Write the pdb file out as a temp file for PyMol FIXME
        pdb = self.database.extractPDBFile(self.pdbref)
        if not pdb:
            pdb = self.database.fetchPDBFileFromWeb(self.pdbref)
        with open("temp.pdb", mode='w') as flines:
            flines.write("\n".join(pdb))

        pymolScript = "\n".join([
            "load temp.pdb, {0}_{1}".format(self.pdbref, i)
            for i in range(len(pymolCommands))
        ])
        pymolScript += "\n"
        pymolScript += "\n".join(pymolCommands)
        pymolScript += """
#formatting
bg_color white
hide all
#show sticks
show cartoon
spectrum b, rainbow,  minimum=0, maximum=1
set opaque_background=0
set antialias = on
set line_smooth = 1
set depth_cue = 1
set specular = 1
set surface_quality = 1
set stick_quality = 15
set sphere_quality = 2
set ray_trace_fog = 0.8
set light = (-0.2,0,-1)

set ray_shadows, 0
set surface_mode, 1
set cartoon_side_chain_helper,on
rebuild
        """

        # If we are doing a single-chain analysis, cut out the other chains
        try:
            chainRef = self.database.extractDocumentGivenId(self.edgelistid)['chainref']
            pymolScript += f"""
select notGivenChain, ! chain {chainRef}
remove notGivenChain
zoom
"""
        except KeyError:
            pass

        # If we are after a png, then generate one
        if outputPng:
            pymolScript += f"""
set ray_trace_mode = 1
png {self.pdbref}.png, width=10cm, dpi=300, ray=1
"""

        with open("temp.pml", mode='w') as flines:
            flines.write(pymolScript)

        if not outputPng:
            subprocess.run(["pymol", "temp.pml"])
        else:
            # Run quietly
            subprocess.run(["pymol", "-c", "temp.pml"])
        os.remove("temp.pml")
        os.remove("temp.pdb")

    def draw(self):
        """Draw the underlying network as a NetworkX graph, colour by community."""
        edgelist = self.database.extractDocumentGivenId(
            self.edgelistid)['data']
        G = nx.Graph()
        for i, j, weight in edgelist:
            G.add_edge(i, j, weight=weight)

        pos = nx.spring_layout(G)
        fig, ax = plt.subplots(figsize=(5, 5))

        # Suppress MPL's complaining, as it's a NetworkX problem.
        warnings.filterwarnings("ignore")
        count = 0
        partition = np.asarray(self.data[0], dtype=int)
        for com in set(partition):  # This will break for a 1d partition
            count = count + 1
            list_nodes = np.where(partition == com + 1)[0]
            nx.draw_networkx_nodes(
                G, pos, list(list_nodes), node_size=20, node_color=Set3_12.mpl_colors[count % 12])

        nx.draw_networkx_edges(G, pos, alpha=0.5)
        ax.set_title("Network for {}".format(self.pdbref))
        plt.show()


def toTree(data):
        """
        Take the partition, and output a tree in which each node is a community, with edge
        weight the community size, connected to its parent node.
        """
        # The first layer is just each top-level community connected to the root node.
        treeDepth = len(data)
        edges = []
        print(data)
        for i in set(data[0]):
            edges.append([0, i])
        print()
        # Now for each sub-row, get the communities belonging to the parent row
        if treeDepth >= 2:
            for i in range(treeDepth - 1):

                column = np.asarray(data[i], dtype=int)
                print(column)
                for community in set(column):
                    # Get the slice of the arrays below correspond to that parent community
                    print(column[column == community])
                    print()
        print()
        print(edges)


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
