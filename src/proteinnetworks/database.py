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
        N : the number of iterations of Infomap run

if doctype == pdbfile:
    data: The PDBfile itself (sans headers) as an array of strings


NB this might be too large for MongoDB to handle (>16MB)
"""

import pymongo
import datetime
import urllib.request
import logging

from pymongo.errors import ConnectionFailure, OperationFailure, DuplicateKeyError
from bson.errors import InvalidId
from bson.objectid import ObjectId

# import sys


class Database:
    """A wrapper around MongoDB."""

    def __init__(self, password="", local=False):
        """Connect to MongoDB, and ensure that it's running."""
        # Server is stored locally; if I can't find it in 1 second its not running.
        if not local:
            if not password:
                password = input("password: ").strip()
            self.client = pymongo.MongoClient(
                "mongodb://writeAccess:" + password +
                "@s7.tcm.phy.private.cam.ac.uk/proteinnetworks",
                serverSelectionTimeoutMS=1000)
            try:
                # The ismaster command is cheap and does not require auth.
                self.client.admin.command('ismaster')
            except ConnectionFailure as err:
                # Propagate the exception back up to whoever called it
                raise IOError("Couldn't connect to the database") from err
            except OperationFailure as err:
                raise IOError("Password incorrect") from err

            self.db = self.client.proteinnetworks
            self.collection = self.db.proteinnetworks
        else:
            self.collection = LocalCollection()

    def extractEdgelist(self,
                        pdbref,
                        edgelisttype,
                        hydrogenstatus,
                        scaling,
                        chainref=None):
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
            "scaling": scaling
        }
        if chainref is not None:
            query["chainref"] = chainref
        else:
            # Needs to explicitly look for the record without a chainref field
            query["chainref"] = {"$exists": False}

        self.validateEdgelist(query, excludeData=True)
        cursor = self.collection.find(query)
        numresults = cursor.count()
        if not numresults:
            return
        elif numresults == 1:
            return cursor[0]
        else:
            raise IOError("More than one edgelist found matching the query")

    def depositEdgelist(self,
                        pdbref,
                        edgelisttype,
                        hydrogenstatus,
                        scaling,
                        edges,
                        chainref=None):
        """
        Deposit edgelist into the database.

        Check that the given edgelist isn't already in the database,
        then deposit and return the _id.
        """
        edgelist = {
            "pdbref": pdbref,
            "doctype": "edgelist",
            "edgelisttype": edgelisttype,
            "hydrogenstatus": hydrogenstatus,
            "scaling": scaling
        }
        if chainref is not None:
            edgelist["chainref"] = chainref
            cursor = self.collection.find(edgelist)

        else:
            # Explictly pass a "doesn't have a chainref field" to the query
            temp = edgelist.copy()
            temp['chainref'] = {"$exists": False}
            cursor = self.collection.find(edgelist)

        numresults = cursor.count()
        if numresults:
            raise IOError(
                "Edgelist already exists in the database! Something has gone terribly wrong!"
            )
        else:
            edgelist["date"] = datetime.datetime.utcnow()
            edgelist["data"] = edges

            self.validateEdgelist(edgelist)

            logging.info("adding edgelist to database...")
            result = self.collection.insert_one(edgelist)
            return result.inserted_id

    def extractPDBFile(self, pdbref):
        """
        Validate the PDB reference and attempt to extract the PDB file corresponding to the PDB ref.

        return None if the PDB reference cannot be found.
        throw an IOError if the PDB reference is invalid
        """
        if not (type(pdbref) == str and len(pdbref) == 4):
            raise IOError("Malformed PDB reference:", pdbref)

        query = {
            "pdbref": pdbref,
            "doctype": "pdbfile",
        }
        cursor = self.collection.find(query)
        numresults = cursor.count()
        if numresults > 1:
            raise IOError("More than one PDB file found")
        elif numresults == 1:
            return cursor[0]['data']

    def fetchPDBFileFromWeb(self, pdbref):
        """
        Pull the PDB file from the web, deposit, and return.

        Fetched from the RCSB, and the headers stripped to save disk space
        """
        # Validate the pdbref
        if not (type(pdbref) == str and len(pdbref) == 4):
            raise IOError("Malformed PDB reference:", pdbref)

        # Pull the PDB file from the Central Repo
        headers = [
            "HEADER", "OBSLTE", "TITLE", "SPLT", "CAVEAT", "COMPND", "SOURCE",
            "KEYWDS", "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT",
            "SPRSDE", "JRNL", "REMARKS", "REMARK"
        ]
        url = "http://www.rcsb.org/pdb/files/{}.pdb".format(pdbref)
        data = urllib.request.urlopen(url).readlines()
        pdbfile = []
        for line in data:
            line = line.decode().strip()
            # Strip out the headers.
            for header in headers:
                if line.startswith(header):
                    break
            else:
                pdbfile.append(line)

        document = {
            "pdbref": pdbref,
            "doctype": "pdbfile",
            "data": pdbfile,
        }
        logging.info("adding PDB file to database...")
        try:
            self.collection.insert_one(document)
        except DuplicateKeyError as err:
            raise IOError("PDB file already in the database") from err
        return pdbfile

    def extractPartition(self, pdbref, edgelistid, detectionmethod, r, N):
        """
        Validate the parameter set and attempt to extract the partition.

        Return None if the parameters are valid, but the partition isn't found.
        """
        # Check that the edgelistid maps to a database entry
        try:
            numberOfEdgelists = self.collection.find({
                "_id":
                ObjectId(edgelistid),
                "doctype":
                "edgelist"
            }).count()
        except InvalidId as err:
            raise IOError("edgelistid not valid:", edgelistid) from err

        if numberOfEdgelists:
            query = {
                "pdbref": pdbref,
                "doctype": "partition",
                "edgelistid": edgelistid,
                "detectionmethod": detectionmethod
            }
            if r != -1:
                query['r'] = r
            if N != -1:
                query['N'] = N

            self.validatePartition(query, excludeData=True)

            cursor = self.collection.find(query)
            numresults = cursor.count()
            if not numresults:
                return
            elif numresults == 1:
                return cursor[0]
            else:
                raise IOError("More than one partition found")
        else:
            logging.error("No edgelist found with the given id")

    def extractDocumentGivenId(self, documentid):
        """Return a document given an id. Return None if not found."""
        try:
            return self.collection.find_one({"_id": ObjectId(documentid)})
        except InvalidId as err:
            raise IOError("Invalid ID") from err

    def depositPartition(self, pdbref, edgelistid, detectionmethod, r, N,
                         data):
        """
        Deposit partition into the database.

        Check that the given partition isn't already in the database,
        then deposit and return the _id.
        """
        partition = {
            "pdbref": pdbref,
            "doctype": "partition",
            "detectionmethod": detectionmethod,
            "edgelistid": edgelistid,
        }

        if r != -1:
            partition['r'] = r
        if N != -1:
            partition['N'] = N

        cursor = self.collection.find(partition)
        numresults = cursor.count()
        if numresults:
            raise IOError(
                "Partition already exists in the database! Something has gone terribly wrong!"
            )
        else:
            partition["date"] = datetime.datetime.utcnow()
            partition["data"] = data

            self.validatePartition(partition)

            logging.info("adding partition to database...")

            result = self.collection.insert_one(partition)
            return result.inserted_id

    def validatePartition(self, partition, excludeData=False):
        """
        Test that the proposed partition has the correct arguments.

        If the partition is invalid, throw an IOError. Otherwise, return None.
        """
        # Validate the pdbref field (string of length 4):
        if type(partition['pdbref']) != str and len(partition['pdbref']) != 4:
            raise IOError("PDB reference malformed: ", partition['pdbref'])
        # Validate the edgelistId (valid ObjectId and refers to edgelist with same pdbref)
        if type(partition['edgelistid']) != ObjectId:
            raise IOError("edgelistid must be of type ObjectId")
        else:
            doc = self.extractDocumentGivenId(partition['edgelistid'])
            if not doc:
                raise TypeError("edgelist referenced does not exist")
            if doc['pdbref'] != partition['pdbref']:
                raise IOError(
                    "edgelist referenced corresponds to different PDB file.")
        # Validate the detection method (AFG or Infomap)
        if partition['detectionmethod'] != "AFG" and partition['detectionmethod'] != "Infomap":
            raise IOError(
                "Only AFG and Infomap are permitted as detection methods.")
        # Validate the r and N values (r float, N integer, must be one or the other)
        if 'r' in partition.keys():
            if "N" in partition.keys():
                raise IOError("r and N cannot both be community parameters")
            if type(partition['r']) != float:
                raise IOError("r must be of type 'float'")
        elif 'N' not in partition.keys():
            raise IOError("Must have one of 'r' and 'N'")
        elif type(partition['N']) != int:
            raise IOError("N must be of type 'int'")
        # Validate the partition itself (array of integers without gaps.)
        if not excludeData:
            if type(partition['data']) != list:
                raise IOError('partition must be a list')
            else:
                # partition['data'] should either be a list of lists or a list.
                for item in partition['data']:
                    if type(item) != int and type(item) != list:
                        raise IOError(
                            'partition must be a list of ints (perhaps nested, not',
                            type(item))
                # if a 1D list
                if not any(isinstance(i, list) for i in partition['data']):
                    numcoms = len(set(partition['data']))
                    try:
                        for i in range(numcoms):
                            # Check there is at least one item in the list for all coms.
                            partition['data'].index(i + 1)
                    except ValueError as err:
                        raise IOError(
                            'partition invalid: gaps found in labelling'
                        ) from err
                # if a nested list
                elif all(isinstance(i, list) for i in partition['data']):
                    for column in partition['data']:
                        numcoms = len(set(column))
                        try:
                            for j in range(numcoms):
                                # Check there is at least one item in the list for all coms.
                                column.index(j + 1)
                        except ValueError as err:
                            raise IOError(
                                'partition invalid: gaps found in labelling'
                            ) from err

    def getNumberOfDocuments(self):
        """Return the total number of documents in the collection."""
        return self.collection.count()

    def validateEdgelist(self, edgelist, excludeData=False):
        """
        Test that "edgelist" has correct field: value formatting.

        Throw an IOError if any field is found to be invalid.
        """
        pass
        # Validate PDB reference
        if type(edgelist['pdbref']) != str or len(edgelist['pdbref']) != 4:
            raise IOError("PDB reference malformed: ", edgelist['pdbref'])
        # Validate edgelisttype
        if (edgelist['edgelisttype'] != "atomic" and edgelist['edgelisttype'] != "residue"):
            raise IOError("edgelisttype must be either 'atomic' or 'residue'")
        # Validate hydrogenstatus
        if (edgelist['hydrogenstatus'] != "noH" and
                edgelist['hydrogenstatus'] != "Hatoms" and
                edgelist['hydrogenstatus'] != "Hbonds"):
            raise IOError(
                "hydrogenstatus must be either 'noH', 'Hatoms' or 'Hbonds'")
        # Validate scaling:
        if type(edgelist['scaling']) != float or edgelist['scaling'] < 0.0:
            raise IOError("scaling must be a non-negative float")
        # Validate the chain reference
        # if edgelist.get("chainref") is not None and type(edgelist["chainref"]) != str:
        #     raise IOError("chain reference, if it exists, must be a string.")
        """
        Validate edges:
        Array correct shape, correct indexing, no self-loops.
        """
        if not excludeData:
            edges = edgelist['data']
            smallestNode = 1000
            for edge in edges:

                if len(edge) != 3:
                    raise IOError("Edgelist has the wrong shape")
                if type(edge[0]) != int or type(edge[1]) != int:
                    raise IOError("Nodes should be integers")
                minNode = min(edge[0], edge[1])
                if minNode < smallestNode:
                    smallestNode = minNode
                if edge[0] == edge[1]:
                    raise IOError("No self-loops permitted")
            if smallestNode != 1:
                raise IOError("Node labelling should start at 1")

    def extractMappings(self, pdbref, mappingtype):
        """Find all documents corresponding to mappings (e.g. to PFAM) for a given pdb ref."""
        # Check that the edgelistid maps to a database entry
        query = {
            "pdbref": pdbref,
            "doctype": "mapping",
            "mappingtype": mappingtype,
        }

        cursor = self.collection.find(query)
        return cursor

    def extractSuperNetwork(self, pdbref, partitionid, level):
        """
        Attempt to extract the supernetwork.

        Return None if the parameters are valid, but the partition isn't found.
        """
        # Check that the edgelistid maps to a database entry
        try:
            numberOfPartitions = self.collection.find({
                "_id":
                ObjectId(partitionid),
                "doctype":
                "partition"
            }).count()
        except InvalidId as err:
            raise IOError("edgelistid not valid:", partitionid) from err

        if numberOfPartitions:
            query = {
                "pdbref": pdbref,
                "doctype": "supernetwork",
                "partitionid": partitionid
            }
            if level is not None:
                query['level'] = level
            cursor = self.collection.find(query)
            numresults = cursor.count()
            if not numresults:
                return
            elif numresults == 1:
                return cursor[0]
            else:
                raise IOError("More than one partition found")
        else:
            logging.error("No partition found with the given id")

    def depositSuperNetwork(self, pdbref, partitionid, level, data):
        """
        Deposit supernetwork into the database.

        Check that the given partition isn't already in the database,
        then deposit and return the _id.
        """
        supernetwork = {
            "pdbref": pdbref,
            "doctype": "supernetwork",
            "partitionid": partitionid,
            "level": level
        }

        cursor = self.collection.find(supernetwork)
        numresults = cursor.count()
        if numresults:
            raise IOError(
                "Partition already exists in the database! Something has gone terribly wrong!"
            )
        else:
            supernetwork['data'] = data

            logging.info("adding supernetwork to database...")

            result = self.collection.insert_one(supernetwork)
            return result.inserted_id

    def extractAllSuperNetworks(self, pdbref=None):
        """Extract all supernetworks, except the one specified by pdbref."""
        query = {
            "pdbref": {
                "$ne": pdbref
            },
            "doctype": "supernetwork",
        }
        cursor = self.collection.find(query)
        return cursor


class LocalCollection:
    """
    A very limited in-memory "database" to be used if no MongoDB instance can be found.

    Stores records as a list of dicts, manipulated with the following methods:
    - collection.find():
        given a set of parameters (including wildcards such as $exists etc) as a dict,
        return a list of the dicts in the database matching this description.

    - collection.find_one():
        as above, but only return one. Used in this code to pull things by their ObjectId.

    - collection.insert_one():
        given a dict, add an ObjectId, push the record, return the id.

    - count():
        return the number of records in the db.
    """

    def __init__(self):
        """Initialise the empty list of dicts."""
        self.storageList = []

    def find(self, query):
        """
        Return a 'cursor' which behaves like a generator with a count method.
        """

        class Cursor(list):
            """Extend the list class with a count method that does the same thing as len()."""

            def count(self):
                return len(self)

        subset = []
        for record in self.storageList:
            for key, value in query.items():
                if key not in record:
                    break
                if type(value) == dict and "$exists" in value:
                    exists = value["$exists"]
                    # match if "exists" is False and key isn't in the record
                    # or if "exists" is True and key is in the record
                    match = (exists and key in record) or (not exists and
                                                           (key not in record))
                    if not match:
                        break
                else:
                    if record[key] != value:
                        break
            else:
                subset.append(record)

        results = Cursor(subset)
        return results

    def find_one(self, query):
        """
        Return a single record.

        As we expect the local DB to be small, this can just be find() with a length check.
        """
        results = self.find(query)
        assert len(results) < 2
        if len(results) == 1:
            return results[0]

    def insert_one(self, record):
        """
        Push a dictionary to the "database", adding a BSON ObjectId, and return a Result
        (with an inserted_id attribute).

        If the record is already there, throw an IOError
        """

        class Result:
            """A container for the inserted_id, necessary to match the pymongo collection."""

            def __init__(self, id):
                self.inserted_id = id

        record["_id"] = ObjectId()
        self.storageList.append(record)
        result = Result(record["_id"])
        return result

    def count(self):
        return len(self.storageList)
