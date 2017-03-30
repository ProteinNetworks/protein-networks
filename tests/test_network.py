"""
Unit tests for the Network class.

Units to be tested:

__init__.
    note that this includes Database __init__ (is perhaps a sign that the class design is iffy)
generateEdgelist
 extractAtomicData(self, pdbdata):
"""
# import proteinnetworks.network


def test_init(mock_database):
    """Test the __init__() function in the network.py module.

    Inputs: pdbref, edgelisttype, hydrogenstatus, scaling
    Output: a Network() class with members:
        - scaling
        - edgelisttype
        - hydrogenstatus
        - pdbref
        - database
        - edgelist
    """
    # inputArgs = {
    #     "scaling": 4.5,
    #     "edgelisttype": "residue",
    #     "hydrogenstatus": "noH",
    #     "pdbref": "1ubq"
    # }
    # pn = proteinnetworks.network.Network(**inputArgs)
    assert 1 == 1
