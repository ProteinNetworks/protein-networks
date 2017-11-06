"""
Unit tests for the Database class.

Functions:
__init__ x
extractEdgelist
depositEdgelist x
extractPDBFile x
fetchPDBFileFromWeb
extractPartition x
extractDocumentGivenId x
depositPartition x
validatePartition
getNumberOfDocuments x
validateEdgelist
extractMappings
extractSuperNetwork
depositSuperNetwork
extractAllSuperNetworks
"""

import proteinnetworks.database
from bson.objectid import ObjectId
import pytest


def test_database_initialise(mock_database):
    """
    Meta-test, to check that pytest and Travis are behaving.

    As the MongoClient interface is mocked for these unit tests,
    the __init__ function doesn't really do anything.
    """
    db = proteinnetworks.database.Database(password="bla")
    assert db


def test_database_getnumberofdocuments(mock_database):
    """
    Test the getNumberOfDocuments function.

    Again, this is a simple wrapper on a mocked class, so hard
    to test.
    """
    db = proteinnetworks.database.Database(password="bla")
    assert db.getNumberOfDocuments() == 5


"""
Test the depositPartition function.

Function inputs: pdbref, edgelistid, detectionmethod, r, N, data.
Function outputs: _id

What is the expected behaviour:
    - When all the arguments are correct? x
    - When the arguments are correct but the partition already
        exists in the database? x
    - When the arguments don't make sense? (Lots of options for this)
"""


def test_database_depositpartition_success(mock_database):
    """Test the depositPartition function when all arguments are correct."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        "pdbref": "1ubq",
        "edgelistid": ObjectId("58dbe03fef677d54224a01da"),
        "detectionmethod": "Infomap",
        "r": -1,
        "N": 5,
        "data": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    }
    resultId = db.depositPartition(**depositionArgs)
    assert type(resultId) == ObjectId


def test_database_depositpartition_doc_already_present(mock_database):
    """Test what happens if you try to add a document that's already there."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': '1ubq',
        'data': [[
            3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 3, 3, 3, 3, 1, 1, 1, 1, 1,
            1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 2, 4, 4,
            4, 4, 4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 4,
            4, 4, 2, 2, 2, 2, 2, 2, 2, 2
        ], [
            37, 38, 35, 34, 68, 70, 71, 74, 75, 76, 73, 72, 69, 42, 36, 45, 39,
            12, 10, 17, 5, 8, 2, 9, 62, 60, 57, 65, 61, 58, 59, 66, 63, 64, 67,
            23, 28, 26, 29, 21, 19, 22, 47, 49, 46, 55, 56, 54, 53, 50, 14, 11,
            18, 7, 6, 3, 15, 13, 1, 16, 4, 44, 43, 41, 40, 52, 48, 51, 20, 27,
            25, 24, 30, 31, 32, 33
        ]],
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_malformed_pdb(mock_database):
    """Assert that depositPartition handles bad PDB references correctly."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "blablabla",
        'data': [1, 2, 3],
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_edgelistid_invalid(mock_database):
    """Assert that an malformed ObjectId is rejected correctly."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "1ubq",
        'data': [1, 2, 3],
        'N': 10,
        'edgelistid': 'bla',
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_edgelistid_incorrect(mock_database):
    """
    Test valid, but incorrect ObjectIDs.

    Assert that if the ObjectId doesn't correspond to an edgelist
    with the same PDB reference as the partition, the system throws an exception.
    """
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "1ubq",
        'data': [1, 2, 3],
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_detectionmethod_invalid(mock_database):
    """Assert that if "detectionmethod" isn't "AFG" or "Infomap" then an exception is raised."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 2, 3],
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'bla',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_r_invalid(mock_database):
    """Test that invalid r values (aka not a float) are treated properly."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 2, 3],
        'N': -1,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': "ten"
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_N_invalid(mock_database):
    """Test that invalid N values (aka not an integer) are treated properly."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 2, 3],
        'N': "five",
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_r_and_N_not_both_missing(mock_database):
    """Test that if r and N are simultaneously -1 an error is raised."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 2, 3],
        'N': -1,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_r_and_N_not_both_present(mock_database):
    """Test that if r and N are simultaneously not -1 an error is raised."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 2, 3],
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': 1.5
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_data_invalid(mock_database):
    """Test that if the partition isn't an array of ints, an error is raised."""
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': "bla",
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


def test_database_depositpartition_data_incorrect(mock_database):
    """
    Test that the system rejects incorrect partitions.

    Partitions should have at least one value in 1,...,M where M is the
    number of unique communities.
    """
    db = proteinnetworks.database.Database(password="bla")

    depositionArgs = {
        'pdbref': "2vcr",
        'data': [1, 1, 3, 4],
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.depositPartition(**depositionArgs)


"""
Test the extractDocumentGivenId function
input - an ObjectId
output - a document

Tests:
        - Does it return the right thing given a correct id?
        - Does it fail sensibly given an incorrect id? (i.e one that doesn't exist in the database)
        - Does it fail sensibly given a super-incorrect id?

Again, a lightweight wrapper around a mocked method, so can't test too hard.
"""


def test_database_extractdocumentgivenid_normal(mock_database):
    """Assert that the database returns the right document given an id."""
    db = proteinnetworks.database.Database(password="bla")
    doc = db.extractDocumentGivenId("58dbe03fef677d54224a01d9")
    assert doc


def test_database_extractdocumentgivenid_missing(mock_database):
    """Assert that the database returns None if the document to be found is missing."""
    db = proteinnetworks.database.Database(password="bla")
    doc = db.extractDocumentGivenId("58dbe03fef677d54224a01d7")
    assert not doc


def test_database_extractdocumentgivenid_malformed(mock_database):
    """Assert that the database returns an IOError if the id is malformed."""
    db = proteinnetworks.database.Database(password="bla")
    with pytest.raises(IOError):
        db.extractDocumentGivenId("bla")


"""
Test the extractPartition function.

Function inputs: pdbref, edgelistid, detectionmethod, r, N.
Function outputs: partition (as dict)

What is the expected behaviour:
    - When all the arguments are correct?
    - When the arguments are correct but the partition doesn't exist
    - When the arguments don't make sense? (Lots of options for this)
"""


def test_database_extractpartition_success(mock_database):
    """Assert that a partition we know to exist can be successfully pulled."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    doc = db.extractPartition(**extractionArgs)
    assert doc['data'] == [[
        3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4,
        4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 4, 4, 4, 2, 2, 2, 2,
        2, 2, 2, 2
    ], [
        37, 38, 35, 34, 68, 70, 71, 74, 75, 76, 73, 72, 69, 42, 36, 45, 39, 12,
        10, 17, 5, 8, 2, 9, 62, 60, 57, 65, 61, 58, 59, 66, 63, 64, 67, 23, 28,
        26, 29, 21, 19, 22, 47, 49, 46, 55, 56, 54, 53, 50, 14, 11, 18, 7, 6,
        3, 15, 13, 1, 16, 4, 44, 43, 41, 40, 52, 48, 51, 20, 27, 25, 24, 30,
        31, 32, 33
    ]]


def test_database_extractpartition_missing(mock_database):
    """Test that a partition we know not to exist (but with correct formatting) returns None."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224b12da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    doc = db.extractPartition(**extractionArgs)
    assert not doc


def test_database_extractpartition_malformed_pdb(mock_database):
    """Assert that depositPartition handles bad PDB references correctly."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': 'blablabla',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.extractPartition(**extractionArgs)


def test_database_extractpartition_edgelistid_invalid(mock_database):
    """Assert that an malformed ObjectId is rejected correctly."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': 'blablabla',
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.extractPartition(**extractionArgs)


def test_database_extractpartition_edgelistid_incorrect(mock_database):
    """
    Test valid, but incorrect ObjectIDs.

    Assert that if the ObjectId doesn't correspond to an edgelist
    with the same PDB reference as the partition, the system throws an exception.
    """
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.extractPartition(**extractionArgs)


def test_database_extractpartition_detectionmethod_invalid(mock_database):
    """Assert that if "detectionmethod" isn't "AFG" or "Infomap" then an exception is raised."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'blab',
        'r': -1
    }
    with pytest.raises(IOError):
        db.extractPartition(**extractionArgs)


def test_database_extractpartition_r_invalid(mock_database):
    """Test that invalid r values (aka not a float) are treated properly."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': b'5'
    }
    with pytest.raises(IOError):
        db.extractPartition(**extractionArgs)


def test_database_extractpartition_N_invalid(mock_database):
    """Test that invalid N values (aka not an integer) are treated properly."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '1ubq',
        'N': 10.0,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.extractPartition(**extractionArgs)


def test_database_extractpartition_r_and_N_not_both_missing(mock_database):
    """Test that if r and N are simultaneously -1 an error is raised."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '1ubq',
        'N': -1,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        db.extractPartition(**extractionArgs)


def test_database_extractpartition_r_and_N_not_both_present(mock_database):
    """Test that if r and N are simultaneously not -1 an error is raised."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': 10
    }
    with pytest.raises(IOError):
        db.extractPartition(**extractionArgs)


"""
Test the extractPDBFile function.

Function inputs: pdbref
Function outputs: the PDB file, with headers stripped.
"""


def test_database_extractpdbfile_success(mock_database):
    """Assert that a correct PDB reference returns the correct PDB file, correctly formatted."""
    db = proteinnetworks.database.Database(password="bla")
    pdbref = "1ubq"
    pdbfile = db.extractPDBFile(pdbref)
    assert pdbfile == [
        'DBREF  1UBQ A    1    76  UNP    P62988   UBIQ_HUMAN       1     76',
        'SEQRES   1 A   76  MET GLN ILE PHE VAL LYS THR LEU THR GLY LYS THR ILE',
        'SEQRES   2 A   76  THR LEU GLU VAL GLU PRO SER ASP THR ILE GLU ASN VAL',
        'SEQRES   3 A   76  LYS ALA LYS ILE GLN ASP LYS GLU GLY ILE PRO PRO ASP',
        'SEQRES   4 A   76  GLN GLN ARG LEU ILE PHE ALA GLY LYS GLN LEU GLU ASP',
        'SEQRES   5 A   76  GLY ARG THR LEU SER ASP TYR ASN ILE GLN LYS GLU SER',
        'SEQRES   6 A   76  THR LEU HIS LEU VAL LEU ARG LEU ARG GLY GLY',
        'FORMUL   2  HOH   *58(H2 O)',
        'HELIX    1  H1 ILE A   23  GLU A   34  1                                  12',
        'HELIX    2  H2 LEU A   56  TYR A   59  5                                   4',
        'SHEET    1 BET 5 GLY A  10  VAL A  17  0',
        'SHEET    2 BET 5 MET A   1  THR A   7 -1',
        'SHEET    3 BET 5 GLU A  64  ARG A  72  1',
        'SHEET    4 BET 5 GLN A  40  PHE A  45 -1',
        'SHEET    5 BET 5 LYS A  48  LEU A  50 -1',
        'CRYST1   50.840   42.770   28.950  90.00  90.00  90.00 P 21 21 21    4',
        'ORIGX1      1.000000  0.000000  0.000000        0.00000',
        'ORIGX2      0.000000  1.000000  0.000000        0.00000',
        'ORIGX3      0.000000  0.000000  1.000000        0.00000',
        'SCALE1      0.019670  0.000000  0.000000        0.00000',
        'SCALE2      0.000000  0.023381  0.000000        0.00000',
        'SCALE3      0.000000  0.000000  0.034542        0.00000',
        'ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N',
        'ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C',
        'ATOM      3  C   MET A   1      26.913  26.639   3.531  1.00  9.62           C',
        'ATOM      4  O   MET A   1      27.886  26.463   4.263  1.00  9.62           O',
        'ATOM      5  CB  MET A   1      25.112  24.880   3.649  1.00 13.77           C',
        'ATOM      6  CG  MET A   1      25.353  24.860   5.134  1.00 16.29           C',
        'ATOM      7  SD  MET A   1      23.930  23.959   5.904  1.00 17.17           S',
        'ATOM      8  CE  MET A   1      24.447  23.984   7.620  1.00 16.11           C',
        'ATOM      9  N   GLN A   2      26.335  27.770   3.258  1.00  9.27           N',
        'ATOM     10  CA  GLN A   2      26.850  29.021   3.898  1.00  9.07           C',
        'ATOM     11  C   GLN A   2      26.100  29.253   5.202  1.00  8.72           C',
        'ATOM     12  O   GLN A   2      24.865  29.024   5.330  1.00  8.22           O',
        'ATOM     13  CB  GLN A   2      26.733  30.148   2.905  1.00 14.46           C',
        'ATOM     14  CG  GLN A   2      26.882  31.546   3.409  1.00 17.01           C',
        'ATOM     15  CD  GLN A   2      26.786  32.562   2.270  1.00 20.10           C',
        'ATOM     16  OE1 GLN A   2      27.783  33.160   1.870  1.00 21.89           O',
        'ATOM     17  NE2 GLN A   2      25.562  32.733   1.806  1.00 19.49           N',
        'ATOM     18  N   ILE A   3      26.849  29.656   6.217  1.00  5.87           N',
        'ATOM    220  O   LYS A  29      38.020  29.772  10.382  1.00  6.87           O',
        'ATOM    221  CB  LYS A  29      36.193  27.058   9.911  1.00 10.28           C',
        'ATOM    222  CG  LYS A  29      36.153  25.620   9.409  1.00 14.94           C',
        'ATOM    223  CD  LYS A  29      34.758  25.280   8.900  1.00 19.69           C',
        'ATOM    224  CE  LYS A  29      34.793  24.264   7.767  1.00 22.63           C',
        'ATOM    594  N   GLY A  75      41.165  35.531  31.898  0.25 36.31           N',
        'ATOM    595  CA  GLY A  75      41.845  36.550  32.686  0.25 36.07           C',
        'ATOM    596  C   GLY A  75      41.251  37.941  32.588  0.25 36.16           C',
        'ATOM    597  O   GLY A  75      41.102  38.523  31.500  0.25 36.26           O',
        'ATOM    598  N   GLY A  76      40.946  38.472  33.757  0.25 36.05           N',
        'ATOM    599  CA  GLY A  76      40.373  39.813  33.944  0.25 36.19           C',
        'ATOM    600  C   GLY A  76      40.031  39.992  35.432  0.25 36.20           C',
        'ATOM    601  O   GLY A  76      38.933  40.525  35.687  0.25 36.13           O',
        'ATOM    602  OXT GLY A  76      40.862  39.575  36.251  0.25 36.27           O',
        'TER     603      GLY A  76',
        'HETATM  604  O   HOH A  77      45.747  30.081  19.708  1.00 12.43           O',
        'HETATM  605  O   HOH A  78      19.168  31.868  17.050  1.00 12.65           O',
        'HETATM  606  O   HOH A  79      32.010  38.387  19.636  1.00 12.83           O',
        'HETATM  659  O   HOH A 132      38.363  30.369   5.579  0.49 35.45           O',
        'HETATM  660  O   HOH A 133      27.841  46.062  17.589  0.81 32.15           O',
        'HETATM  661  O   HOH A 134      37.667  43.421  17.000  0.50 33.32           O',
        'MASTER      274    0    0    2    5    0    0    6  660    1    0    6',
        'END'
    ]


def test_database_extractpdbfile_missing(mock_database):
    """Assert that if the PDB reference is valid but not found, nothing is returned."""
    db = proteinnetworks.database.Database(password="bla")
    pdbref = "1abc"
    pdbfile = db.extractPDBFile(pdbref)
    assert not pdbfile


def test_database_extractpdbfile_error(mock_database):
    """Assert that if the PDB reference is invalid, an IOError is thrown."""
    db = proteinnetworks.database.Database(password="bla")
    pdbref = 5
    with pytest.raises(IOError):
        db.extractPDBFile(pdbref)


"""
Test the fetchPDBFileFromWeb function.

Function inputs: pdbref
Function outputs: the PDB file, with headers stripped.
Also should deposit the PDB file.
"""


def test_database_fetchpdbfilefromweb_success(mock_database, mock_urlopen):
    """Assert that a valid PDB reference can be fetched, stripped, deposited and returned."""
    db = proteinnetworks.database.Database(password="bla")
    pdbref = "3rty"
    pdbfile = db.fetchPDBFileFromWeb(pdbref)
    assert len(pdbfile) == 11


def test_database_fetchpdbfilefromweb_error(mock_database, mock_urlopen):
    """Assert that a malformed PDB reference will cause an IOError."""
    db = proteinnetworks.database.Database(password="bla")
    pdbref = 5
    with pytest.raises(IOError):
        db.fetchPDBFileFromWeb(pdbref)


def test_database_fetchpdbfilefromweb_pdbfilealreadyindatabase(mock_database, mock_urlopen):
    """
    Try to fetch a PDB file that's already in the database.

    Assert than an IOError is thrown.
    """
    db = proteinnetworks.database.Database(password="bla")
    pdbref = "1ubq"
    with pytest.raises(IOError):
        db.fetchPDBFileFromWeb(pdbref)


"""
Test the depositEdgelist function.

Function inputs: pdbref, edgelisttype, hydrogenstatus, scaling, edges.
Function outputs: _id

What is the expected behaviour:
    - When all the arguments are correct?
    - When the arguments are correct but the partition already
        exists in the database?
    - When the arguments don't make sense? (Lots of options for this)
"""


def test_database_depositedgelist_normal(mock_database):
    """Test that if all arguments are correct, the edgelist is deposited succesfully."""
    db = proteinnetworks.database.Database(password="bla")
    depositionArgs = {
        'pdbref': '2vc5',
        'edges':
        [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]],
        'edgelisttype': 'residue',
        'hydrogenstatus': 'noH',
        'scaling': 4.5
    }
    resultId = db.depositEdgelist(**depositionArgs)
    assert type(resultId) == ObjectId


def test_database_depositedgelist_already_present(mock_database):
    """Test that if the edgelist is already present, an exception is thrown."""
    db = proteinnetworks.database.Database(password="bla")
    depositionArgs = {
        'pdbref': '2vcr',
        'edges':
        [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]],
        'edgelisttype': 'residue',
        'hydrogenstatus': 'noH',
        'scaling': 4.5
    }
    with pytest.raises(IOError):
        db.depositEdgelist(**depositionArgs)


def test_database_depositedgelist_pdbref_invalid(mock_database):
    """Test that a malformed PDB reference is rejected."""
    db = proteinnetworks.database.Database(password="bla")
    depositionArgs = {
        'pdbref': 'blablabla',
        'edges':
        [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]],
        'edgelisttype': 'residue',
        'hydrogenstatus': 'noH',
        'scaling': 4.5
    }
    with pytest.raises(IOError):
        db.depositEdgelist(**depositionArgs)


def test_database_depositedgelist_edgelisttype_invalid(mock_database):
    """Test that if the edgelisttype isn't "atomic" or "residue" it is rejected."""
    db = proteinnetworks.database.Database(password="bla")
    depositionArgs = {
        'pdbref': '2vc5',
        'edges':
        [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]],
        'edgelisttype': 'blablabla',
        'hydrogenstatus': 'noH',
        'scaling': 4.5
    }
    with pytest.raises(IOError):
        db.depositEdgelist(**depositionArgs)


def test_database_depositedgelist_hydrogenstatus_invalid(mock_database):
    """Test that if the hydrogenstatus isn't "noH", "Hatoms" or "Hbonds" it is rejected."""
    db = proteinnetworks.database.Database(password="bla")
    depositionArgs = {
        'pdbref': '2vc5',
        'edges':
        [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]],
        'edgelisttype': 'residue',
        'hydrogenstatus': 'blablabla',
        'scaling': 4.5
    }
    with pytest.raises(IOError):
        db.depositEdgelist(**depositionArgs)


def test_database_depositedgelist_scaling_invalid_negative(mock_database):
    """Assert that scaling values must be 0.0."""
    db = proteinnetworks.database.Database(password="bla")
    depositionArgs = {
        'pdbref': '2vc5',
        'edges':
        [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]],
        'edgelisttype': 'residue',
        'hydrogenstatus': 'noH',
        'scaling': -5.0
    }
    with pytest.raises(IOError):
        db.depositEdgelist(**depositionArgs)


def test_database_depositedgelist_scaling_invalid_notnumber(mock_database):
    """Assert that scaling values must be floats."""
    db = proteinnetworks.database.Database(password="bla")
    depositionArgs = {
        'pdbref': '2vc5',
        'edges':
        [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]],
        'edgelisttype': 'residue',
        'hydrogenstatus': 'noH',
        'scaling': "five"
    }
    with pytest.raises(IOError):
        db.depositEdgelist(**depositionArgs)


def test_database_depositedgelist_edges_not_array(mock_database):
    """Assert that the edgelist is rejected if it is not a _ x 3 array of numbers."""
    db = proteinnetworks.database.Database(password="bla")
    depositionArgs = {
        'pdbref': '2vc5',
        'edges':
        [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70, 102], [5, 3, 23]],
        'edgelisttype': 'residue',
        'hydrogenstatus': 'noH',
        'scaling': 4.5
    }
    with pytest.raises(IOError):
        db.depositEdgelist(**depositionArgs)


def test_database_depositedgelist_edges_not_correctly_indexed(mock_database):
    """Assert that the edge labelling starts from 1."""
    db = proteinnetworks.database.Database(password="bla")
    depositionArgs = {
        'pdbref': '2vc5',
        'edges':
        [[2, 6, 44], [3, 6, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]],
        'edgelisttype': 'residue',
        'hydrogenstatus': 'noH',
        'scaling': 4.5
    }
    with pytest.raises(IOError):
        db.depositEdgelist(**depositionArgs)


def test_database_depositedgelist_edges_no_self_loops(mock_database):
    """Assert no self-loops."""
    db = proteinnetworks.database.Database(password="bla")
    depositionArgs = {
        'pdbref': '2vc5',
        'edges':
        [[1, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]],
        'edgelisttype': 'residue',
        'hydrogenstatus': 'noH',
        'scaling': 4.5
    }
    with pytest.raises(IOError):
        db.depositEdgelist(**depositionArgs)


"""
Test the extractEdgelist function.

Function inputs: pdbref, edgelisttype, hydrogenstatus, scaling
Function outputs: edgelist (as dict)

What is the expected behaviour:
    - When all the arguments are correct?
    - When the arguments are correct but the partition doesn't exist
    - When the arguments don't make sense? (Lots of options for this)
"""


def test_database_extractedgelist_success(mock_database):
    """Ensure successful extraction if all arguments correct."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '2vcr',
        'edgelisttype': 'residue',
        'hydrogenstatus': 'noH',
        'scaling': 4.5
    }
    doc = db.extractEdgelist(**extractionArgs)
    assert doc['data'] == [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]]


def test_database_extractedgelist_edgelist_missing(mock_database):
    """Ensure None returned if edgelist missing."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': '2vcr',
        'edgelisttype': 'residue',
        'hydrogenstatus': 'Hatoms',
        'scaling': 4.5
    }
    doc = db.extractEdgelist(**extractionArgs)
    assert not doc


def test_database_extractedgelist_pdbref_invalid(mock_database):
    """Ensure exception thrown if pdbref malformed."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': 5,
        'edgelisttype': 'residue',
        'hydrogenstatus': 'Hatoms',
        'scaling': 4.5
    }
    with pytest.raises(IOError):
        db.extractEdgelist(**extractionArgs)


def test_database_extractedgelist_edgelisttype_invalid(mock_database):
    """Ensure exception thrown if edgelisttype malformed."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': 5,
        'edgelisttype': 'residues',
        'hydrogenstatus': 'Hatoms',
        'scaling': 4.5
    }
    with pytest.raises(IOError):
        db.extractEdgelist(**extractionArgs)


def test_database_extractedgelist_hydrogenstatus_invalid(mock_database):
    """Ensure exception thrown if hydrogenstatus not from accepted options."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': 5,
        'edgelisttype': 'residue',
        'hydrogenstatus': [],
        'scaling': 4.5
    }
    with pytest.raises(IOError):
        db.extractEdgelist(**extractionArgs)


def test_database_extractedgelist_scaling_invalid(mock_database):
    """Ensure exception thrown if scaling invalid."""
    db = proteinnetworks.database.Database(password="bla")
    extractionArgs = {
        'pdbref': 5,
        'edgelisttype': 'residue',
        'hydrogenstatus': 'Hatoms',
        'scaling': -4.5
    }
    with pytest.raises(IOError):
        db.extractEdgelist(**extractionArgs)
