"""Database interface methods (perhaps should be a class).

MongoDB database format:
database: proteinnetworks
collection: proteinnetworks
document has fields:
    _id: bla
    pdbref: bla
    doctype: (edgelist | partition | pdbfile ) currently

if doctype == edgelist:
    edgelisttype: (atomic | residue)
    hydrogenstatus: (noH | Hatoms | Hbonds) (i.e. preprocessed with nothing, HAAD, or Stride)
    scaling: the value of the cutoff used
    data: the actual edgelist

if doctype == partition:
    edgelistid: The _id of the edgelist used in generating the partition
    detectionmethod: (AFG | Infomap) as it stands

    if detectionmethod == AFG:
        r : The AFG parameter used in generating the partition
        data: The 1D array giving the partition

    if detectionmethod == Infomap:
        data: A 2D numpy array giving the data

if doctype == pdbfile:
    data: The PDBfile itself


NB this might be too large for MongoDB to handle (>16MB)
"""

import pymongo
from pymongo.errors import ConnectionFailure
import datetime
import urllib.request

# import sys


class Database:
    """A wrapper around MongoDB."""

    def __init__(self):
        """Connect to MongoDB, and ensure that it's running."""
        # Server is stored locally; if I can't find it in 1 second its not running.
        self.client = pymongo.MongoClient(serverSelectionTimeoutMS=1000)
        try:
            # The ismaster command is cheap and does not require auth.
            self.client.admin.command('ismaster')
        except ConnectionFailure:
            # Propagate the exception back up to whoever called it
            raise IOError

        self.db = self.client.proteinnetworks
        self.collection = self.db.proteinnetworks

    def extractEdgelist(self, pdbref, edgelisttype, hydrogenstatus, scaling):
        """
        Attempt to extract the edgelist matching the given parameter set.

        Return None if the edgelist cannot be found.

        There should never be two documents with the same parameter combination (for now)
        """
        query = {
            "pdbref": pdbref,
            "doctype": "edgelist",
            "edgelisttype": edgelisttype,
            "hydrogenstatus": hydrogenstatus,
            "scaling": scaling,
        }
        cursor = self.collection.find(query)
        numresults = cursor.count()
        if not numresults:
            return
        elif numresults == 1:
            return cursor[0]
        else:
            raise IOError  # TODO Custom exception

    def depositEdgelist(self,
                        pdbref,
                        edgelisttype,
                        hydrogenstatus,
                        scaling,
                        edges,
                        dryrun=True):
        """Check that the given edgelist isn't already in the database, then deposit."""
        edgelist = {
            "pdbref": pdbref,
            "doctype": "edgelist",
            "edgelisttype": edgelisttype,
            "hydrogenstatus": hydrogenstatus,
            "scaling": scaling,
        }
        cursor = self.collection.find(edgelist)
        numresults = cursor.count()
        if numresults:
            print(
                "Edgelist already exists in the database! Something has gone terribly wrong!"
            )
        else:
            edgelist["date"] = datetime.datetime.utcnow()
            edgelist["data"] = edges
            print("adding to database...")
            if not dryrun:
                self.collection.insert_one(edgelist)
            else:
                print("not really")

    def extractPDBFile(self, pdbref):
        """
        Pull the PDB file corresponding to the given PDB ref, if it exists.

        If it doesn't exist, pull it from the web, strip out the headers, and
        add it to the database.
        """
        headers = [
            "HEADER",
            "OBSLTE",
            "TITLE",
            "SPLT",
            "CAVEAT",
            "COMPND",
            "SOURCE",
            "KEYWDS",
            "EXPDTA",
            "NUMMDL",
            "MDLTYP",
            "AUTHOR",
            "REVDAT",
            "SPRSDE",
            "JRNL",
            "REMARKS",
        ]
        query = {
            "pdbref": pdbref,
            "doctype": "pdbfile",
        }
        cursor = self.collection.find(query)
        numresults = cursor.count()
        if numresults > 1:
            raise IOError
        elif numresults == 1:
            return cursor[0]['data']
        else:
            # Pull the PDB file from the Central Repo
            url = "http://www.rcsb.org/pdb/files/{}.pdb".format(pdbref)

            print("Fetching {}.pdb...".format(pdbref))
            data = urllib.request.urlopen(url).readlines()
            print(data)
            # Strip out the headers.