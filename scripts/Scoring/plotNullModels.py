"""
Plot a given partition, and its null models.
"""
from getAllJaccardZscore import generateNullModel
import numpy as np
import proteinnetworks
from bson.objectid import ObjectId
import matplotlib.pyplot as plt
from palettable.colorbrewer.qualitative import Set3_12


with open("../pdbswithstorededgelists.dat") as flines:
    pdbs = [line.strip() for line in flines]

pdb = np.random.choice(pdbs).strip()

db = proteinnetworks.database.Database()
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

# This is the partition to be tested
testpartition = partition.data[maxI]


# Generate an arbitrary number of null models, and plot them
numModels = 6
nullJaccard = []
for i in range(numModels):
    nullmodel = generateNullModel(partition.data[maxI])
    nullJaccard.append([proteinnetworks.insight.getModifiedJaccard(pfamDomains, nullmodel), nullmodel])

numPlots = 1 + numModels
fig, axes = plt.subplots(nrows=numPlots, figsize=(5, 8), sharex=True)
data = np.vstack(x[1] for x in nullJaccard)

axes[0].imshow(
    np.vstack(
        2 * [testpartition]),  # vstack otherwise imshow complains
    aspect=50,
    cmap=Set3_12.mpl_colormap)
axes[0].xaxis.set_ticks_position('bottom')
axes[0].yaxis.set_visible(False)



for i, ax in enumerate(axes[1:]):
    ax.imshow(
        np.vstack(
            2 * [data[i, :]]),  # vstack otherwise imshow complains
        aspect=50,
        cmap="viridis")
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_visible(False)

plt.suptitle(pdb)
plt.tight_layout()
plt.subplots_adjust(top=0.95)
plt.xlabel('Residue number')
plt.savefig("NullModels/{}.png".format(pdb), dpi=300)
plt.show()