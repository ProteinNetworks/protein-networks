"""Plot the modified Jaccard for the ~8000 Infomap results with scaling 4.5 on residue networks."""

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

sns.set_context("poster")

with open("allJaccards.dat") as flines:
    # get the average TMScore for each line
    contactJaccards = [
        float(x.strip().split(" ")[-1]) for x in flines
        if len(x.strip().split(" ")) > 1
    ]

contactJaccards = np.asarray(contactJaccards)

fig, ax = plt.subplots()
sns.distplot(contactJaccards, kde=False)  # label="Isomorphs")
ax.set_xlim([0, 1])
ax.set_ylim(bottom=0)
plt.xlabel('Modified Jaccard')
plt.ylabel('Frequency Density', labelpad=20)

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(30)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(25)

plt.savefig("contactJaccards.png", dpi=300)
plt.show()
