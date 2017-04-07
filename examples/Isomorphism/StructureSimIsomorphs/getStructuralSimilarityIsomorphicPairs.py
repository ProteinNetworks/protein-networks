"""Get the TM scores for a sample of the base data, and the isomorphic data."""
import subprocess
import re
import itertools
import proteinnetworks
import numpy as np

# Read the original data
with open("../isomorphicProteins.dat") as flines:
    isomorphicClasses = [
        line.strip().split(" ") for line in flines
        if len(line.strip().split(" ")) != 1
    ]

# Convert to pairs
isomorphicPairs = []
for isomorphicClass in isomorphicClasses:
    for pair in itertools.combinations(isomorphicClass, 2):
        isomorphicPairs.append(pair)


sampleSize = 2000
# print(isomorphicPairs)
# Choose 500 of the isomorphic pairs
indices = np.random.randint(len(isomorphicPairs), size=sampleSize)
sampledIsomorphicPairs = []
for index in indices:
    sampledIsomorphicPairs.append(list(isomorphicPairs[index]))

# Choose 500 pairs from the base data
db = proteinnetworks.database.Database()

baserate = db.collection.find({"doctype": "supernetwork"})
baserate = [x['pdbref'] for x in baserate]
basePairs = []
for i in range(sampleSize):
    pairs = np.random.choice(baserate, size=2)
    basePairs.append(pairs.tolist())

with open("isomorphicPairs.dat", mode='w') as flines:
    flines.write(repr(sampledIsomorphicPairs))

with open("basePairs.dat", mode='w') as flines:
    flines.write(repr(basePairs))

sampledIsomorphicPairsStructureSim = []
for pair in sampledIsomorphicPairs:
    pdb1, pdb2 = pair
    # Get the PDB files from the database, retrieve the TM-score

    pdbfile1 = db.extractPDBFile(pdb1)
    pdbfile2 = db.extractPDBFile(pdb2)
    with open("temp1.pdb", mode='w') as flines:
        flines.write("\n".join(line for line in pdbfile1))
    with open("temp2.pdb", mode='w') as flines:
        flines.write("\n".join(line for line in pdbfile2))

    # Calculate the TM score
    output = subprocess.run([
        "/home/will/MainProject/External Code/TMalign/TMalign", "temp1.pdb",
        "temp2.pdb", "-a"
    ], stdout=subprocess.PIPE).stdout
    TMScore = re.findall(b"(?<=TM-score\= )\d\.\d*", output)
    if TMScore:
        TMScore1, TMScore2 = [float(x) for x in TMScore]
        # print(TMScore1, TMScore2)
        sampledIsomorphicPairsStructureSim.append(
            [pdb1, pdb2, str(TMScore1), str(TMScore2)])

with open("isomorphicPairsWithStructureSim.dat", mode='w') as flines:
    flines.write("\n".join(" ".join(line)
                           for line in sampledIsomorphicPairsStructureSim))

basePairsStructureSim = []
for pair in basePairs:
    pdb1, pdb2 = pair
    # Get the PDB files from the database, retrieve the TM-score

    pdbfile1 = db.extractPDBFile(pdb1)
    pdbfile2 = db.extractPDBFile(pdb2)
    with open("temp1.pdb", mode='w') as flines:
        flines.write("\n".join(line for line in pdbfile1))
    with open("temp2.pdb", mode='w') as flines:
        flines.write("\n".join(line for line in pdbfile2))

    # Calculate the TM score
    output = subprocess.run([
        "/home/will/MainProject/External Code/TMalign/TMalign", "temp1.pdb",
        "temp2.pdb", "-a"
    ], stdout=subprocess.PIPE).stdout

    TMScore = re.findall(b"(?<=TM-score\= )\d\.\d*", output)
    if TMScore:
        TMScore1, TMScore2 = [float(x) for x in TMScore]
        # print(TMScore1, TMScore2)
        basePairsStructureSim.append(
            [pdb1, pdb2, str(TMScore1), str(TMScore2)])

with open("baseRateWithStructureSim.dat", mode='w') as flines:
    flines.write("\n".join(" ".join(line)
                           for line in basePairsStructureSim))
