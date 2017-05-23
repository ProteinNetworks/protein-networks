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
    prevI = -1
    numBoundaries = -1  # Start from -1 to avoid counting the start as a boundary
    for i in testPartition:
        if i != prevI:
            numBoundaries += 1
        prevI = i

    # Place the boundaries arbitrarily, and assign each segment to a randomly chosen community
    newBoundaries = []
    for i in range(numBoundaries):
        randint = np.random.randint(low=1, high=len(testPartition))
        while randint in newBoundaries:
            randint = np.random.randint(len(testPartition))
        newBoundaries.append(randint)
    newBoundaries.sort()
    newBoundaries = [0] + newBoundaries + [len(testPartition) - 1]
    # Fill the gaps with communities in the range 1... numComs.

    # We want to randomly assign, but ensure the same number of communities.
    # So assign minimum number of communities, then add to reach quota.
    newCommunities = [i + 1 for i in range(numCommunities)]
    while len(newCommunities) != len(newBoundaries) - 1:
        newCommunities.append(np.random.randint(numCommunities) + 1)

    # Shuffle until no two numbers are together
    np.random.shuffle(newCommunities)
    valid = False
    while not valid:
        valid = True
        np.random.shuffle(newCommunities)
        for i in range(len(newCommunities) - 1):
            if newCommunities[i] == newCommunities[i + 1]:
                valid = False
        

    nullModel = np.zeros(len(testPartition))
    for i in range(len(newBoundaries) - 1):
        nullModel[newBoundaries[i]:newBoundaries[i + 1] + 1] = newCommunities[i]

    # Check the null model
    # Get the number of boundaries
    prevI = -1
    nullModelNumBoundaries = -1  # Start from -1 to avoid counting the start as a boundary
    for i in nullModel:
        if i != prevI:
            nullModelNumBoundaries += 1
        prevI = i
    assert set(nullModel) == set(testPartition)
    assert len(nullModel) == len(testPartition)   
    if nullModelNumBoundaries != numBoundaries:
        np.set_printoptions(threshold=np.nan)
        print(newCommunities)
        print(newBoundaries)
        print(nullModel)
        print()
        print(testPartition)
        raise NotImplementedError
    return nullModel

if __name__ == "__main__":

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
                                                        
                print(nullJaccard)
                aosidjasoidj
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
