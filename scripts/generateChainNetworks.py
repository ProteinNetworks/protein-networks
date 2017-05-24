"""
generateChainNetworks.py

For all proteins with partitions in the DB, create networks corresponding
to each chain.
"""
import proteinnetworks
from bson.objectid import ObjectId


db= proteinnetworks.database.Database()

with open("./pdbswithstoredpartitions.dat") as flines:
    pdbs = [x.strip() for x in flines]

size = len(pdbs)
for i, pdb in enumerate(pdbs):
    print(f"{i} of {size} completed")
    pdbdata = db.extractPDBFile(pdb)
            
    if not pdbdata:
        pdbdata = db.fetchPDBFileFromWeb(pdb)
    pdbChains = set(x[21] for x in pdbdata if x[:4] == "ATOM")
    for chainref in pdbChains:
        inputArgs = {"scaling": 4.0,
                "edgelisttype": "residue",
                "hydrogenstatus": "noH",
                "pdbref": pdb,
                "chainref": chainref,
                "database": db}
        network = proteinnetworks.network.Network(**inputArgs)
        partitionArgs = {"pdbref": pdb,
                         "edgelistid": ObjectId(network.edgelistid),
                         "detectionmethod": "Infomap",
                        "N": 1000,
                        "database": db}
        partition = proteinnetworks.partition.Partition(**partitionArgs)
    
