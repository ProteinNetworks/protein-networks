"""
Generate Jaccard indices for the generated partitions. Output to a file.

Then get the z-score for each entry:

z-score = ( J - mu / sigma ) i.e. the number of standard deviations away from the mean.
Gives the significance.
Plot as a histogram.
"""
import proteinnetworks
import numpy as np
from bson.objectid import ObjectId

with open("../pdbswithstorededgelists.dat") as flines:
    pdbs = [line.strip() for line in flines]

db = proteinnetworks.database.Database()

totalpdbs = len(pdbs)
jaccards = []
try:
    for i, pdb in enumerate(pdbs):
        try:
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

            # Get the Jaccard
            pfamDomains = np.asarray(partition.getPFAMDomainArray(), dtype=int)
            maxJaccard = -1
            maxI = -1

            for i, col in enumerate(partition.data):
                jaccard = proteinnetworks.insight.getModifiedJaccard(
                    pfamDomains, np.asarray(
                        col, dtype=int))
                if jaccard > maxJaccard:
                    maxJaccard = jaccard
                    maxI = i
            if maxI != -1:
                jaccards.append([pdb, maxI, maxJaccard])
        except FileNotFoundError as e:
            print(str(e))
        except ValueError as e:
            print(str(e))
except KeyboardInterrupt:
    with open("allJaccards.dat", mode='w') as flines:
        flines.write("\n".join(" ".join(map(str, x)) for x in jaccards))

with open("allJaccards.dat", mode='w') as flines:
    flines.write("\n".join(" ".join(map(str, x)) for x in jaccards))
