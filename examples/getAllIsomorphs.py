"""Get all isomorphic supernetworks."""

import proteinnetworks
from tqdm import tqdm
import os

db = proteinnetworks.database.Database()

supernetworks = db.collection.find({"doctype": "supernetwork"})
supernetworks = list(supernetworks)
if os.path.exists("isomorphicProteins.dat"):
    with open("isomorphicProteins.dat") as flines:
        isomorphismClasses = [line.strip().split(" ") for line in flines]
else:
    isomorphismClasses = [[]]
try:
    for supernetwork in tqdm(supernetworks):
        match = False
        # Check whether the protein has already been assigned a class
        for isomorphismClass in isomorphismClasses:
            if supernetwork['pdbref'] in isomorphismClass:
                match = True
                break
        if match:
            continue
        # If not already in a class, generate the SuperNetwork
        supernetwork = proteinnetworks.insight.SuperNetwork.fromPartitionId(
            supernetwork['partitionid'], db)

        isomorphs = supernetwork.getIsomorphs()
        isomorphs.append(supernetwork.pdbref)
        isomorphismClasses.append(isomorphs)
    with open("isomorphicProteins.dat", mode='w') as flines:
        flines.write("\n".join(" ".join(x) for x in isomorphismClasses if x))
except KeyboardInterrupt:
    with open("isomorphicProteins.dat", mode='w') as flines:
        flines.write("\n".join(" ".join(x) for x in isomorphismClasses if x))
