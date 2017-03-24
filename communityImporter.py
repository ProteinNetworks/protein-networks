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
        N: the infomap N i.e. number of repeated runs (if known)

if doctype == pdbfile:
    data: The PDBfile itself


NB this might be too large for MongoDB to handle (>16MB)
"""

import pymongo
import glob
import os
import datetime
from pythonModules import treeFileToNumpyArray

if __name__ == "__main__":
    client = pymongo.MongoClient()
    db = client.proteinnetworks
    collection = db.proteinnetworks

    # The contact networks are stored in /home/will/MainProject/Topology/contactnetworks/????/,
    #  where ???? is the PDB reference

    treefiles = glob.glob(
        "/home/will/MainProject/Topology/contactnetworks/????/????.?.?.tree")
    # edgelisttype = "residue"
    doctype = "partition"
    detectionmethod = "Infomap"
    N = -1  # I don't know what Infomap N I used to generate this data
    unassignedTreeFiles = []
    for treefile in treefiles:
        pdbref = os.path.basename(treefile)[:4]
        print(pdbref)
        scaling = float(os.path.basename(treefile)[5:-5])
        # Get the objectID for the edgelist used in generating the partition
        cursor = collection.find({
            "doctype": "edgelist",
            "pdbref": pdbref,
            "scaling": scaling
        })
        if cursor.count() != 1:
            unassignedTreeFiles.append(treefile)
            continue
        edgelistid = cursor[0]["_id"]
        data = treeFileToNumpyArray(treefile)
        edgelistDocument = {
            "pdbref": pdbref,
            "doctype": doctype,
            "edgelistid": edgelistid,
            "detectionmethod": detectionmethod,
            "data": data.tolist(),
            "N": N,
            "date": datetime.datetime.utcnow(),
        }
        result = collection.insert_one(edgelistDocument)
        print(result.inserted_id)
