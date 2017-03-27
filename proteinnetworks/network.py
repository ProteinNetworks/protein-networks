"""Stores functionality related to the generation and analysis of edgelists."""

from .database import Database
# from database import Database
import sys


class Network:
    """Holds a edgelist and its parameters, and offers network inspection methods."""

    def __init__(self, pdbref, edgelisttype, hydrogenstatus, scaling):
        """
        Initialise the edgelist with a given parameter set.

        Edgelist specified by:
            - scaling i.e cutoff
            - atomic / residue
            - status of hydrogen atoms
            - PDB reference
        """
        self.scaling = scaling
        self.edgelisttype = edgelisttype
        self.hydrogenstatus = hydrogenstatus
        self.pdbref = pdbref

        # Try to connect to the database
        try:
            database = Database()
            print("successfully connected")
        except IOError:
            print("Couldn't connect to server")
            sys.exit()
        # Attempt to extract the edgelist matching the given params
        doc = database.getEdgelist(pdbref, edgelisttype, hydrogenstatus, scaling)
        if doc:
            self.edgelist = doc['data']
            print("edgelist found")
        else:
            print("no edgelist fitting those parameters found: generating")