"""
Creates the MongoDB database, from all the data I have scattered around my filesystem.

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
import glob
import os
import datetime

if __name__ == "__main__":
    client = pymongo.MongoClient()
    db = client.proteinnetworks
    collection = db.proteinnetworks

    # The contact networks are stored in /home/will/MainProject/Topology/contactnetworks/????/,
    #  where ???? is the PDB reference

    edgelists = glob.glob(
        "/home/will/MainProject/Topology/contactnetworks/????/????.?.?.dat")
    edgelisttype = "residue"
    doctype = "edgelist"
    hydrogenstatus = "noH"
    for edgelist in edgelists:
        pdbref = os.path.basename(edgelist)[:4]
        scaling = float(os.path.basename(edgelist)[5:-4])
        lastmodified = datetime.datetime.utcnow()
        assert scaling in [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        with open(edgelist) as flines:
            # edgelist is a list of lines of the form (node1 node2 weight)
            edges = []
            for line in flines:
                node1, node2, weight = line.strip().split(" ")
                edges.append([node1, node2, weight])
        # Create the document to be stored
        edgelistDocument = {
            "pdbref": pdbref,
            "doctype": doctype,
            "edgelistype": edgelisttype,
            "hydrogenstatus": hydrogenstatus,
            "scaling": scaling,
            "date": datetime.datetime.utcnow(),
            "data": edges
        }
        # Add the data to the "proteinnetworks" collection
        result = collection.insert_one(edgelistDocument)
        print(result)
