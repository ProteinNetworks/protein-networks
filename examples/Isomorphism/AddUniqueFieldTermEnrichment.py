# Add 'isUnique' to the GO terms JSON

import json
with open("./GOtermsCorrected.json") as flines:
    enrichment = json.load(flines)

GOterms = []
for isoClass in enrichment:
    print("pdbs:")
    print(" ".join(isoClass['isoClass']))
    print("GO terms with p <0.01:")
    enrichedTerms = [x for x in isoClass['GOterms'] if float(x['pCorr']) < 0.01]
    # Sort by enrichment
    enrichedTerms = sorted(enrichedTerms, key=lambda k: float(k['pCorr']))
    for term in enrichedTerms:
        print(term['GOlabel'], term['pCorr'])
    GOterms.append([[term['GOlabel'], term['pCorr']] for term in enrichedTerms])
    print()
    print()

GOtermsSet = []
for isoClass in GOterms:
    GOtermsSet.append(set(x[0] for x in isoClass))


for i, isoClass in enumerate(GOtermsSet):

    otherSets = set.union(*[elem for j, elem in enumerate(GOtermsSet) if i != j])

    uniques = isoClass - otherSets
    for j in enrichment[i]['GOterms']:
        if j['GOlabel'] in uniques and float(j['pCorr']) < 0.01:
            j['isUnique'] = True
        else:
            j['isUnique'] = False

with open("./GOtermsCorrectedUnique.json", mode='w') as flines:
    json.dump(enrichment, flines)
