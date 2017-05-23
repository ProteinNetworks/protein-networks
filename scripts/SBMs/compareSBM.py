"""Compare the output of SBM with Infomap."""

import glob
import os
from bson.objectid import ObjectId
import numpy as np
import matplotlib.pyplot as plt
from palettable.colorbrewer.qualitative import Set3_12
import proteinnetworks


def getInfomapData(pdbRef):
    """Given a PDB reference, return the list of lists from the database."""
    global db
    inputArgs = {
        "scaling": 4.0,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": pdbRef,
        "database": db
    }
    proteinNetwork = proteinnetworks.network.Network(**inputArgs)
    partitionArgs = {
        "pdbref": pdbRef,
        "edgelistid": ObjectId(proteinNetwork.edgelistid),
        "detectionmethod": "Infomap",
        "N": 100,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    return partition.data


def plotStripeDiagramSideBySide(infomapData, SBMdata, title):
    """Plot the partition as a set of "stripe" plots."""
    # Convert the list (or nested list) to a numpy array, for use with imshow.
    stripes1 = np.asarray(infomapData, dtype=int)
    stripes2 = np.asarray(SBMdata, dtype=int)

    numInfomapArrays = np.shape(stripes1)[0]
    numSBMArrays = np.shape(stripes2)[0]
    while numInfomapArrays < numSBMArrays:

        # Pad with ones until both arrays are the same shape
        stripes1 = np.vstack([stripes1, [1] * np.shape(stripes1)[1]])
        numInfomapArrays = np.shape(stripes1)[0]

    while numSBMArrays < numInfomapArrays:
        # Pad with ones until both arrays are the same shape
        stripes2 = np.vstack([stripes2, [1] * np.shape(stripes2)[1]])
        numSBMArrays = np.shape(stripes2)[0]

    assert numInfomapArrays == numSBMArrays
    numPlots = numInfomapArrays

    fig, axes = plt.subplots(
        nrows=numPlots, ncols=2, figsize=(8, 8), sharex=True)
    try:
        stripes3 = np.stack([stripes1, stripes2], axis=2)
    except ValueError:
        print(np.shape(stripes1), np.shape(stripes2))
    assert np.all(stripes3[:, :, 0] == stripes1)
    for i in range(numInfomapArrays):
        axes[i, 0].imshow(
            np.vstack(
                2 * [stripes3[i, :, 0]]),  # vstack otherwise imshow complains
            aspect=100,
            cmap=Set3_12.mpl_colormap,
            interpolation="nearest")
        axes[i, 0].xaxis.set_ticks_position('bottom')
        axes[i, 0].yaxis.set_visible(False)
    for i in range(numSBMArrays):
        axes[i, 1].imshow(
            np.vstack(
                2 * [stripes3[i, :, 1]]),  # vstack otherwise imshow complains
            aspect=100,
            cmap=Set3_12.mpl_colormap,
            interpolation="nearest")
        axes[i, 1].xaxis.set_ticks_position('bottom')
        axes[i, 1].yaxis.set_visible(False)
    axes[0, 0].set_title("Infomap")
    axes[0, 1].set_title("SBM")
    plt.suptitle(title)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.xlabel('Residue number')
    plt.savefig("Data/{}.png".format(pdbRef), dpi=300)
    plt.clf()

SBMfiles = glob.glob("Data/*.sbm")
db = proteinnetworks.database.Database(password="8moMlmkRnZurXsl")

for SBMfile in SBMfiles:
    with open(SBMfile) as flines:
        # The + 1 is to maintain consistency with Infomap's labelling.
        SBMdata = [[int(x) + 1 for x in line.strip().split(" ")]
                   for line in flines]

        # Get the corresponding Infomap file
        pdbRef = os.path.basename(SBMfile)[:4]
        try:
            Infomapdata = getInfomapData(pdbRef)
        except FileNotFoundError:
            continue
        plotStripeDiagramSideBySide(Infomapdata, SBMdata, pdbRef)
