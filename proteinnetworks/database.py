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
import sys


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
            print("Server not available:")
            print("Log file stored in /var/log/mongod/mongod.log")
            sys.exit(1)

        self.db = self.client.proteinnetworks
        self.collection = self.db.proteinnetworks
