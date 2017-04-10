"""Get all isomorphic supernetworks."""

import proteinnetworks
from tqdm import tqdm
import os
import numpy as np

db = proteinnetworks.database.Database()

supernetworks = db.collection.find({"doctype": "supernetwork"})
supernetworks = list(supernetworks)
supernetworks = np.random.choice(np.asarray(supernetworks), size=1000, replace=False)
if os.path.exists("weaklyIsomorphicProteins.dat"):
    with open("weaklyIsomorphicProteins.dat") as flines:
        isomorphicPairs = [line.strip().split(" ") for line in flines]
else:
    isomorphicPairs = []
try:
    for supernetwork in supernetworks:
        match = False
        # Check whether the protein has already been assigned a class
        for isomorphicPair in isomorphicPairs:
            if supernetwork['pdbref'] in isomorphicPair:
                match = True
                break
        if match:
            continue
        # If not already in a class, generate the SuperNetwork
        supernetwork = proteinnetworks.insight.SuperNetwork.fromPartitionId(
            supernetwork['partitionid'], db)

        # [[pdb1, pdb2, simscore], [pdb1, pdb3, simscore]]
        isomorphs = supernetwork.getWeakIsomorphs(supernetworks)
        isomorphicPairs += isomorphs
    with open("weaklyIsomorphicProteins.dat", mode='w') as flines:
        flines.write("\n".join(" ".join(map(str, x)) for x in isomorphicPairs if x))
except KeyboardInterrupt:
    with open("weaklyIsomorphicProteins.dat", mode='w') as flines:
        flines.write("\n".join(" ".join(map(str, x)) for x in isomorphicPairs if x))

        # flines.write("\n".join(" ".join(map(str, x)) for x in isomorphicPairs if x))
