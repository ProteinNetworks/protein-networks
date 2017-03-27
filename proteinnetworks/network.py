"""Stores functionality related to the generation and analysis of edgelists."""

import sys
import numpy as np
import math
from .database import Database
from .atomicradii import atomicRadii


class Network:
    """Holds a edgelist and its parameters, and offers network inspection methods."""

    def __init__(self, pdbref, edgelisttype, hydrogenstatus, scaling):
        """
        Initialise the edgelist with a given parameter set.

        Edgelist specified by:
            - scaling i.e cutoff
            - atomic / residue
            - status of hydrogen atoms
            - PDB reference
        """
        self.scaling = scaling
        self.edgelisttype = edgelisttype
        self.hydrogenstatus = hydrogenstatus
        self.pdbref = pdbref

        # Try to connect to the database
        try:
            database = Database()
            print("successfully connected")
        except IOError:
            print("Couldn't connect to server")
            sys.exit()
        # Attempt to extract the edgelist matching the given params
        doc = database.extractEdgelist(pdbref, edgelisttype, hydrogenstatus,
                                       scaling)
        if doc:
            self.edgelist = doc['data']
            print("edgelist found")
        else:
            print("no edgelist fitting those parameters found: generating")
            edgelist = self.generateEdgelist(pdbref, edgelisttype,
                                             hydrogenstatus, scaling)
            database.depositEdgelist(pdbref, edgelisttype, hydrogenstatus,
                                     scaling, edgelist)

    def generateEdgelist(self, pdbref, edgelisttype, hydrogenstatus, scaling):
        r"""
        Generate the edgelist using the supplied parameters.

        The given pdbref is read in from the database (if it exists) else an
        attempt is made to acquire it from the web. A network is then generated in
        which atoms correspond to nodes, and edges are generated such that:
        $ A_{ij} = = 1 - \frac{|d_{ij}|}{c_{ij}} $

        Where $ c_{ij} = s(r_i + r_j) $, and $ r_i $ is the atomic radius of atom i.
        s is here the "scaling".

        FIXME this only currently works for single-atom elements in the pdb file (due to HAAD)
        FIXME I should integrate HAAD.
        FIXME I should sort out STRIDE

        """
        edges = []
        # if args.hydrogens:
        #     proc = subprocess.run(["haad", filename], stderr= subprocess.DEVNULL)

        # if proc.returncode == -11:
        #     print("what a janky code. SIGSEGV caught.")
        #     sys.exit(1)
        #     # HAAD generates a *.pdb.h file, in which the positions of the hydrogens have
        #     # been added and all non-ATOM files (and the element info) stripped out.
        #     filename = filename + ".h"

        filename = "bla"

        positions, elements, residues = self.extractAtomicData(filename)
        assert len(positions) == len(residues) == len(elements)

        # get matrix of square distances
        distance_squared = np.sum(
            (positions[:, np.newaxis, :] - positions[np.newaxis, :, :])**2,
            axis=-1)

        n = np.shape(positions)[0]

        if edgelisttype == "atomic":
            # Then use atoms as vertices
            for i in range(0, n):
                for j in range(0, i):
                    # distance = getDistance(positions, i, j)
                    cutoff = (atomicRadii[elements[i]] +
                              atomicRadii[elements[j]]) * scaling
                    # cutoff = cutoffs[i,j]
                    if distance_squared[i, j] < cutoff * cutoff:
                        weight = (cutoff - math.sqrt(distance_squared[i, j])
                                  ) / cutoff
                        edges.append([i + 1, j + 1, weight])

        elif edgelisttype == "residue":
            # use the residues as vertices
            edgeList = {}
            for i in range(0, n):
                for j in range(0, i):
                    cutoff = (atomicRadii[elements[i]] +
                              atomicRadii[elements[j]]) * scaling
                    if distance_squared[i, j] < cutoff * cutoff:
                        res1, res2 = residues[i], residues[j]
                        if res1 != res2:
                            if not (res1, res2) in edgeList:
                                edgeList[(res1, res2)] = 1
                            else:
                                edgeList[(res1, res2)] += 1

            # Push to a list and sort by first value
            edges = []
            for row, weight in edgeList.items():
                i, j = row
                edges.append([i, j, weight])
            edges.sort()

    def extractAtomicData(filename):
        """
        Given a path to a filename, extract the atomic data.

        pull all ATOM records, and push the atomic positions, element, and
        residue number to three arrays.
        """
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
                    positions.append(
                        [linelist[30:38], linelist[38:46], linelist[46:54]])
                    # elements.append(linelist[76:78].strip())
                    elements.append(linelist[13].strip())
                    residues.append(residueCounter)
        positions = np.asarray(positions, dtype=float)
        return positions, elements, residues
