"""
getCorrelationOfTMScoreAndMCS.py.

Check whether the network similarity of two proteins correlates with the
TMScore.
"""
import proteinnetworks
import subprocess
import re

with open("weaklyIsomorphicProteinsWithSimScore.dat") as flines:
    weakIsomorphs = [line.strip().split() for line in flines]

db = proteinnetworks.database.Database()

correlation = []

for pdb1, pdb2, networkScore, simScore in weakIsomorphs:

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
        TMScore1, TMScore2 = [x.decode() for x in TMScore]
        # print(TMScore1, TMScore2)
        correlation.append(
            [pdb1, pdb2, networkScore, simScore, TMScore1, TMScore2])

with open("correlationBetweenTMandMCS.dat", mode='w') as flines:
    flines.write("\n".join(" ".join(line)
                           for line in correlation))
