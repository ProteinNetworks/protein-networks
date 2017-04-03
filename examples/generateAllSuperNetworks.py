"""For all proteins with a PFAM mapping, generate their supernetwork."""

import proteinnetworks
from bson.objectid import ObjectId

with open("./pdbswithstorededgelists.dat") as flines:
    pdbs = [line.strip() for line in flines]

db = proteinnetworks.database.Database()

totalpdbs = len(pdbs)

for i, pdb in enumerate(pdbs):
    try:
        print(pdb, i / totalpdbs)
        # Extract/Generate the edgelist
        inputArgs = {
            "scaling": 4.5,
            "edgelisttype": "residue",
            "hydrogenstatus": "noH",
            "pdbref": pdb,
            "database": db
        }
        proteinNetwork = proteinnetworks.network.Network(**inputArgs)
        # Extract/Generate the partition
        partitionArgs = {
            "pdbref": pdb,
            "edgelistid": ObjectId(proteinNetwork.edgelistid),
            "detectionmethod": "Infomap",
            "N": 100,
            "database": db
        }
        partition = proteinnetworks.partition.Partition(**partitionArgs)
        superNetwork = proteinnetworks.insight.SuperNetwork(
            inputpartition=partition)
    except Exception:
        continue
