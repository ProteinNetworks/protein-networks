import proteinnetworks

inputArgs = {
    "scaling": 4.5,
    "edgelisttype": "residue",
    "hydrogenstatus": "noH",
    "pdbref": "16pk"
}
proteinNetwork = proteinnetworks.network.Network(**inputArgs)
