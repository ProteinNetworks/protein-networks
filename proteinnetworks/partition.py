"""Stores functionality related to the generation and analysis of community structures."""

import sys

from .database import Database


class Partition:
    """
    Holds a community structure for a network (i.e the partition) and its parameters.

    Offers partition inspection and visualisation methods.
    """

    def __init__(self, pdbref, edgelistid, detectionmethod, r=-1):
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

        # Try to connect to the database
        try:
            self.database = Database()
            print("successfully connected")
        except IOError:
            print("Couldn't connect to server")
            sys.exit()
        # Attempt to extract the edgelist matching the given params
        doc = self.database.extractPartition(pdbref, edgelistid,
                                             detectionmethod, r)
        if doc:
            self.partition = doc['data']
            print("partition found")
        else:
            print("no partition fitting those parameters found: generating")

            partition = self.generatePartition(pdbref, edgelistid,
                                               detectionmethod, r)
            self.partition = partition
            # self.database.depositPartition(pdbref, edgelistid, detectionmethod,
            #                               r, partition)

    def generatePartition(self, pdbref, edgelistid, detectionmethod, r):
        """Generate a community structure using the parameters supplied."""
        # Get the network.
        edgelist = self.database.extractDocumentGivenId(edgelistid)['data']
        pass
