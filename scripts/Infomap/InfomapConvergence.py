# coding: utf-8

# Check how the number of outer loops (N) changes the conductance
#  (as a proxy for "goodness") of the modules, for a variety of protein shapes and sizes
from bson.objectid import ObjectId
from typing import List
import numpy as np
import proteinnetworks
import json


db = proteinnetworks.database.Database()

# Test Proteins:
# 1fl3 1r76 2ghw 3mym 3rj3 3rpi

for i, pdbRef in enumerate(["1fl3", "1r76", "2ghw", "3mym", "3rj3", "3rpi"]):
    inputArgs = {
        "scaling": 4.0,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": pdbRef,
        "database": db
    }
    proteinNetwork = proteinnetworks.network.Network(**inputArgs)
    convergence = {}
    for N in [1, 10, 100, 1000, 10000]:
        partitionArgs = {
            "pdbref": pdbRef,
            "edgelistid": ObjectId(proteinNetwork.edgelistid),
            "detectionmethod": "Infomap",
            "N": N,
            "database": db
        }
        partition = proteinnetworks.partition.Partition(**partitionArgs)
        conductance = proteinnetworks.insight.getModularity(network=proteinNetwork, partition=partition)
        convergence[N] = conductance

    with open("Qconvergence{}.json".format(i), mode='w') as flines:
        json.dump(convergence, flines, indent=2)
