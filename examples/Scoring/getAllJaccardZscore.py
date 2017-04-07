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


def generateNullModel(testPartition):
    """
    From a given partition, generate a null model.

    Here the null model has the same number of boundaries as the generated partition, but with
    the boundaries arbitrarily placed (and the same number of communities in total.)
    """
    # Get the total number of communities in the partition
    numCommunities = len(set(testPartition))

    # Get the number of boundaries
    prevI = 1
    numBoundaries = -1  # Start from -1 to avoid counting the start as a boundary
    for i in testPartition:
        if i != prevI:
            numBoundaries += 1
        prevI = i

    # Place the boundaries arbitrarily, and assign each segment to a randomly chosen community
    print(len(testPartition))
    newBoundaries = []
    for i in range(numBoundaries):
        newBoundaries.append(np.random.randint(len(testPartition)))
    newBoundaries.sort()
    # Fill the gaps with communities in the range 1... numComs.
    nullModel = np.zeros(len(testPartition))
    print(newBoundaries)
    # for i in range(len(newBoundaries) - 1):
    #     nullModel[newBoundaries[i]:newBoundaries[i + 1] +
    #               1] = i % numCommunities + 1
    # Deal with the start and end bits:

    # partition.plotStripeDiagram()
    np.set_printoptions(threshold=np.nan)
    print(nullModel)
    spifjasopdijf


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

            # Generate an arbitrary number of null models.
            numModels = 100
            nullJaccard = []
            for i in range(numModels):
                nullmodel = generateNullModel(partition.data[maxI])
                nullJaccard.append(
                    proteinnetworks.insight.getModifiedJaccard(pfamDomains,
                                                               nullmodel))

            print(np.mean(nullJaccard))

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
