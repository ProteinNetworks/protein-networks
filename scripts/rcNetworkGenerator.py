"""Parses a Protein Data Bank file into a Numpy Array of the positions.
Generates an edgelist for a residue contact network; i.e. in which each
node is a residue, with edges weighted according to how many links there are
between the residues


"""

import argparse
import numpy as np
import math
import os
from numba import jit


@jit
def getDistance(positionmatrix, i, j):
    """Get the absolute distance in angstroms between residues i and j.

    Arguments:
    positionmatrix - nx3 numpy array with the x,y,z positions of each residue
    i - integer 0 <= i < n
    j - integer 0 <= j < n
    """
    return math.sqrt((positionmatrix[i, 0] - positionmatrix[j, 0])**2 +
                     (positionmatrix[i, 1] - positionmatrix[j, 1])**2 + (positionmatrix[i, 2] - positionmatrix[j, 2])**2)

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="the PDB file to be parsed")
parser.add_argument("-s", "--scaling", type=float, help="the scaling of the cutoff radius (default=1)", default=1)
parser.add_argument("--atoms", help="create graph using atoms", action="store_true")
args = parser.parse_args()
filename = args.filename
scaling = args.scaling
text, anyextension = os.path.splitext(os.path.basename(filename))
output = text + "." + str(scaling) + '.dat'
pos = output + '.pos'
print("Generating edgelist with scaling {0}:".format(scaling))
# For multispecies-cutoff radius; list of the radius of each atom commonly encountered
atomicRadii = dict([('Ac', 2.15), ('Ag', 1.45), ('Al', 1.21), ('Am', 1.8), ('Ar', 1.06),
                    ('As', 1.19), ('At', 1.50), ('Au', 1.36), ('B', 0.84), ('Ba', 2.15),
                    ('Be', 0.96), ('Bh', 1.0), ('Bi', 1.48), ('Bk', 1.0), ('Br', 1.2),
                    ('C', 0.76), ('Ca', 1.76), ('Cd', 1.44), ('Ce', 2.04), ('Cf', 1.0),
                    ('Cl', 1.02), ('Cm', 1.69), ('Co', 1.26), ('Cn', 1.0), ('Cr', 1.39),
                    ('Cs', 2.44), ('Cu', 1.32), ('Db', 1.0), ('Ds', 1.0), ('Er', 1.89),
                    ('Es', 1.0), ('Eu', 1.98), ('F', 0.57), ('Fe', 1.32), ('Fm', 1.0),
                    ('Fr', 2.6), ('Ga', 1.22), ('Gd', 1.96), ('Ge', 1.20), ('H', 0.31),
                    ('He', 0.28), ('Hf', 1.75), ('Hg', 1.32), ('Ho', 1.92), ('Hs', 1.0),
                    ('I', 1.39), ('In', 1.42), ('Ir', 1.41), ('K', 2.03), ('Kr', 1.16),
                    ('La', 2.07), ('Li', 1.28), ('Lr', 1.0), ('Lu', 1.87), ('Md', 1.0),
                    ('Mg', 1.41), ('Mn', 1.39), ('Mo', 1.54), ('Mt', 1.0), ('N', 0.71),
                    ('Na', 1.66), ('Nb', 1.64), ('Nd', 2.01), ('Ne', 0.58), ('Ni', 1.24),
                    ('No', 1.0), ('Np', 1.9), ('O', 0.66), ('Os', 1.44), ('P', 1.07),
                    ('Pa', 2.0), ('Pb', 1.46), ('Pd', 1.39), ('Pm', 1.99), ('Po', 1.40),
                    ('Pr', 2.03), ('Pt', 1.36), ('Pu', 1.87), ('Ra', 2.21), ('Rb', 2.2), ('Re', 1.51),
                    ('Rf', 1.0), ('Rg', 1.0), ('Rh', 1.42), ('Rn', 1.0), ('Ru', 1.46),
                    ('S', 1.05), ('Sb', 1.39), ('Sc', 1.7), ('Se', 1.2), ('Sg', 1.0),
                    ('Si', 1.11), ('Sm', 1.98), ('Sn', 1.39), ('Sr', 1.95), ('Ta', 1.7),
                    ('Tb', 1.94), ('Tc', 1.47), ('Te', 1.38), ('Th', 2.06), ('Ti', 1.6),
                    ('Tl', 1.45), ('Tm', 1.90), ('U', 1.96), ('V', 1.53), ('W', 1.62),
                    ('Xe', 1.40), ('Y', 1.9), ('Yb', 1.87), ('Zn', 1.22), ('Zr', 1.75)])
positions = []
elements = []
residues = []
with open(filename, mode='r') as afile:
    residueCounter = 0
    prevRes = 0
    for line in afile:
        if line.strip() == "ENDMDL":
            break
        linelist = line.rstrip()
        if linelist[0:4] == "ATOM":
            residueNumber = int(linelist[22:26].strip())
            if residueNumber != prevRes:
                residueCounter += 1
            prevRes = residueNumber

            positions.append([linelist[30:38], linelist[38:46], linelist[46:54]])
            elements.append(linelist[76:78].strip())
            residues.append(residueCounter)
positions = np.asarray(positions, dtype=float)
n = np.shape(positions)[0]
# print(residues)
assert len(positions) == len(residues) == len(elements)

"""
for each atom pair, as before:
if the atom pair are within the cutoff distance, then:
    get the residue associated with each atom.
    if the residue numbers are different, then increment the edge weight by 1.

"""
edgeList = {}
for i in range(0, n):
    for j in range(0, i):
        distance = getDistance(positions, i, j)
        cutoff = (atomicRadii[elements[i]] + atomicRadii[elements[j]]) * scaling
        if distance < cutoff:
            res1, res2 = residues[i], residues[j]
            if res1 != res2:
                if not (res1, res2) in edgeList:
                    edgeList[(res1, res2)] = 1
                else:
                    edgeList[(res1, res2)] += 1

# Push to a list and sort by first value
edgeListSorted = []
for row, weight in edgeList.items():
    i, j = row
    edgeListSorted.append([i, j, weight])
edgeListSorted.sort()

with open(output, mode='w') as outputFile:
    for row in edgeListSorted:
        outputFile.write("{} {} {}\n".format(row[0], row[1], row[2]))


print(edgeListSorted)
