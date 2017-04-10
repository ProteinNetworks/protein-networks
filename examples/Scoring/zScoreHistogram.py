"""Plot the z-score associated with the modified Jaccard for the ~8000 residue networks."""

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

sns.set_context("poster")

with open("allZScores.dat") as flines:
    # get the average TMScore for each line
    contactZScores = [x.strip().split(" ")[-1] for x in flines if x.strip().split(" ")[-1] != "nan"]

contactZScores = np.asarray(contactZScores, dtype=float)

fig, ax = plt.subplots()
sns.distplot(contactZScores, kde=True)  # label="Isomorphs")
# ax.set_xlim([0, 1])
# ax.set_ylim(bottom=0)
plt.xlabel('z-score')
plt.ylabel('Frequency Density', labelpad=20)
plt.title("z-score for ~8000 contact networks with scaling 4.5,\n assessed with Infomap")
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(25)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)

plt.savefig("contactZScores.png", dpi=300)
plt.show()
