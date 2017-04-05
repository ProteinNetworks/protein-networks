#!/usr/bin/env python3
"""
Get the Eigenvector Similarity for two weighted edgelists
"""

import proteinnetworks
import networkx as nx


def select_k(spectrum, minimum_energy=0.9):
    running_total = 0.0
    total = sum(spectrum)
    if total == 0.0:
        return len(spectrum)
    for i in range(len(spectrum)):
        running_total += spectrum[i]
        if running_total / total >= minimum_energy:
            return i + 1
    return len(spectrum)


db = proteinnetworks.database.Database()

with open("weaklyIsomorphicProteins2.dat") as flines:
    weakIsomorphs = [line.strip().split(" ") for line in flines]


similarities = []

for pdb1, pdb2, simscore1 in weakIsomorphs:
    edgelist1 = db.collection.find_one({"doctype": "supernetwork", "pdbref": pdb1})['data']
    edgelist2 = db.collection.find_one({"doctype": "supernetwork", "pdbref": pdb2})['data']
    G1 = nx.Graph()
    for i, j, weight in edgelist1:
        G1.add_edge(i, j)
    G2 = nx.Graph()
    for i, j, weight in edgelist2:
        G2.add_edge(i, j)

    laplacian1 = nx.spectrum.laplacian_spectrum(G1, weight=None)
    laplacian2 = nx.spectrum.laplacian_spectrum(G2, weight=None)

    k1 = select_k(laplacian1)
    k2 = select_k(laplacian2)
    k = min(k1, k2)

    similarity = sum((laplacian1[:k] - laplacian2[:k])**2)
    similarities.append([pdb1, pdb2, simscore1, similarity])


with open("simscorecomparison.dat", mode='w') as flines:
    flines.write("\n".join(" ".join(map(str, x)) for x in similarities if x))
