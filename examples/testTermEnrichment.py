"""
Find enriched GO terms within an isomorphism class.

Inputs: - isomorphismClasses.dat, a space-separated list of proteins with isomorphic supernetworks
        - isomorphismClassesGO.json, a list of all GO-terms associated with the PDBs used.

"""
from scipy.stats import hypergeom
import json
import pandas as pd
import numpy as np


def getTermEnrichment(isoClass):
    """
    Take a list of PDBs, output a list of GO terms that occur more often than expected by chance.

    If: n = the number of PDBs is the isoClass
        k = the number of times that a given term occurs
        K = the number of times the term occurs in the base data (aka all proteins
            for which a reduced network was made )
        N = the total number of PDBs in the base data.

    Then the p-value comes from the CDF of the hypergeometric distribution
                                               (binomial, without replacement)

    TODO DEPLETION

    So we want a p-value for each term in the isoClass.
    Steps:
        - Get list of all GO terms in the isoClass
        - Work out K for the base data
        - Work out k for the isoClass
        - Plug into hypergeom.
    """
    n = len(isoClass)
    # basePDBs = glob.glob("../contactnetworks/????/*reduced.dat")
    # print(n)
    global basePDBs
    global goDict
    global df
    global goCounts

    goTermsInIsoClass = {
    }  # Key = GO term. Value = number of times it appears in the isomorphism class
    totalGoTerms = 0
    for pdb in isoClass:
        goTerms = set(goDict[pdb])  # Remove duplicates
        totalGoTerms += len(goTerms)
        for goTerm in goTerms:
            if goTerm in goTermsInIsoClass:
                goTermsInIsoClass[goTerm]['k'] += 1
            else:
                goTermsInIsoClass[goTerm] = {'k': 1}
    # Check your counting is good
    # print(goTermsInIsoClass)
    assert totalGoTerms == sum(x['k'] for x in goTermsInIsoClass.values())

    # drop terms with k = 1
    goTermsInIsoClass = {
        x: y
        for x, y in goTermsInIsoClass.items() if y['k'] != 1
    }

    # Now we have k for each GO Term. Need K
    for goTerm in goTermsInIsoClass:
        goTermsInIsoClass[goTerm]['K'] = int(goCounts[goTerm])
        K = goTermsInIsoClass[goTerm]['K']
        k = goTermsInIsoClass[goTerm]['k']
        if k > 1:
            goTermsInIsoClass[goTerm]['hypergeom_cdf'] = hypergeom.cdf(k, N, n,
                                                                       K)

        else:
            goTermsInIsoClass[goTerm]['hypergeom_cdf'] = -1.0
    # Clean any insignificant terms, take p=0.01 (TODO without Bonferroni currently)
    goTermsInIsoClass = {
        x: y
        for x, y in goTermsInIsoClass.items()
        if y['hypergeom_cdf'] > 0.99 or np.isnan(y['hypergeom_cdf'])
    }
    # goTermsInIsoClass['N'] = N
    # goTermsInIsoClass['n'] = n
    goTermsInIsoClass['PDBsInIsoClass'] = isoClass

    return goTermsInIsoClass


if __name__ == "__main__":

    with open("./isomorphismClasses.dat") as flines:
        isoClasses = [line.strip().split(" ") for line in flines]

    with open("baseData.dat") as flines:
        basePDBs = [line.strip() for line in flines]
    N = len(basePDBs)

    with open("isomorphismClassesGO.json") as dictFile:
        goDict = json.load(dictFile)

    # Pull the data into pandas. Then filter such that only rows with PDB
    # in basePDBs is there (For getting K later).
    df = pd.read_csv(
        "../GO/pdb_chain_go.tsv", sep="\t", header=0, usecols=["PDB", "GO_ID"])
    df = df[df['PDB'].isin(basePDBs)]
    goCounts = df['GO_ID'].value_counts()
    enrichedTerms = []
    for isoClass in isoClasses:
        enrichedTerms.append(getTermEnrichment(isoClass))

    print(json.dumps(enrichedTerms, indent=2, sort_keys=True))
