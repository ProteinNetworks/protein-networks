"""
isomorphicPairsWithStructureSim.dat contains pairs of proteins for which the
module network is isomorphic, with the TMScore for the parent proteins

baseRate.dat gives random pairs from the original dataset, with TMScores

Is isomorphic graphs genuinely exhibit structural similarity, then histograms
should reveal this.
"""

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

with open("isomorphicPairsWithStructureSim.dat") as flines:
    # get the average TMScore for each line
    isomorphicPairs = []
    for line in flines:
        tmscore1, tmscore2 = line.strip().split(" ")[2:]
        isomorphicPairs.append((float(tmscore1) + float(tmscore2)) / 2)

isomorphicPairs = np.asarray(isomorphicPairs)

with open("baseRateWithStructureSim.dat") as flines:
    # get the average TMScore for each line
    baseRate = []
    for line in flines:
        tmscore1, tmscore2 = line.strip().split(" ")[2:]
        baseRate.append((float(tmscore1) + float(tmscore2)) / 2)

baseRate = np.asarray(baseRate)
fig, ax = plt.subplots(nrows=1)
sns.kdeplot(isomorphicPairs, shade=True, label="Isomorphs")
sns.kdeplot(baseRate, shade=True, label="Original Dataset")
plt.legend()
plt.savefig("IsomorphicProteinSimilarity.png", dpi=300)
plt.show()
