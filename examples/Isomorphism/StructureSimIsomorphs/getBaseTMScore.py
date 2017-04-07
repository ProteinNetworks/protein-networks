"""
Get the base rate of TM score, by finding, from all proteins with a .4.5.dat,
the average TM score.
"""

import glob
import numpy as np
import subprocess
import re
# from tqdm import tqdm
proteins = np.asarray(glob.glob("../contactnetworks/????/*.4.5.dat"))

sampleSize = 500
randomProteinSample = np.random.choice(
    proteins, size=sampleSize, replace=False)
TMScores = []
for i, protein in enumerate(randomProteinSample):
    print(i, protein)
    for protein2 in randomProteinSample:
        if protein != protein2:
            # print(protein2)
            pdb1 = protein[:-7] + "pdb"  # the hackiest
            pdb2 = protein2[:-7] + "pdb"
            output = subprocess.run([
                "/home/will/MainProject/External Code/TMalign/TMalign", pdb1,
                pdb2, "-a"
            ], stdout=subprocess.PIPE).stdout
            TMScore = re.findall(b"(?<=TM-score\= )\d\.\d*", output)
            if TMScore:
                TMScore1, TMScore2 = [float(x) for x in TMScore]
                TMScores.append([pdb1[-8:-4], pdb2[-8:-4], TMScore1, TMScore2])

with open("baseRate500.dat", mode='w') as flines:
    flines.write("\n".join(" ".join(str(x) for x in line) for line in TMScores))
