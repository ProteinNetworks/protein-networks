import proteinnetworks.partition
from bson.objectid import ObjectId
import os
import subprocess
import pytest

"""
Unit tests for the partition module.

Units:
Partition:
    __init__
    generatePartition
    plotStripeDiagram
    getPFAMDomainArray
    plotPymolStructure
    draw
treeFileToNestedLists
"""

"""
Tests for Partition.__init__()

inputs: pdbref, edgelistid, detectionmethod, r, N.
outputs: a Partition

Options: - params legit, partition in database
         - params malformed (will be picked up by validatePartition)
         - params legit, partition not in database
"""


def test_partition_init_partition_in_database(mock_database):
    """Test that if the partition is in the DB, it is extracted successfully."""
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    assert type(partition.partitionid) == ObjectId


def test_partition_init_partition_not_in_database(mock_database, mock_subprocess):
    """
    Test that if the partition is not in the DB, it is generated successfully.

    Uses a mocked community detection method.
    """
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 1,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    assert type(partition.partitionid) == ObjectId


# def test_partition_init_partition_in_database_local(mock_database):
#     """Test that if the partition is in the DB, it is extracted successfully."""
#     db = proteinnetworks.database.Database(password="bla")
#     partitionArgs = {
#         'pdbref': '1ubq',
#         'N': 10,
#         'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
#         'detectionmethod': 'Infomap',
#         'r': -1,
#         "database": db
#     }
#     partition = proteinnetworks.partition.Partition(**partitionArgs)
#     assert type(partition.partitionid) == ObjectId


def test_partition_init_partition_no_database_provided(mock_database, mock_subprocess):
    """
    Test that if the partition is not in the DB, it is generated successfully, with
    a tmeporary copy of a database

    Uses a mocked community detection method.
    """

    partitionArgs = {
        'pdbref': '2vc5',
        'N': 1,
        'edgelistid': ObjectId("58dbe03fef677d54224a01da"),
        'detectionmethod': 'Infomap',
        'r': -1
    }
    with pytest.raises(IOError):
        partition = proteinnetworks.partition.Partition(**partitionArgs)
    

"""
Tests for the treeFileToNestedLists function (not a method)

Takes in a path to a treefile, outputs a list of lists.

test: if the treefile doesn't exist we throw a FileNotFoundError
      if the treefile format is unexpected we throw an IOError
      if the treefile is legit we return a list of lists.
"""


def test_partition_treefiletonestedlists_success(mock_database, mock_subprocess):
    """
    Test that a correctly formatted treefile will generate a list of lists.

    Uses the mocked "subprocess.run()" to generate the tree file.
    """
    subprocess.run(["Infomap"])
    listoflists = proteinnetworks.partition.treeFileToNestedLists("temp.tree")
    assert listoflists
    os.remove("temp.tree")


def test_partition_treefiletonestedlists_error(mock_database, mock_subprocess):
    """
    Test that a a FileNotFoundError is thrown if you try to run treeFileToNestedLists
    on something that isn't there
    """
    with pytest.raises(FileNotFoundError):
        proteinnetworks.partition.treeFileToNestedLists("temp2.tree")


"""
tests for getPFAMDomainArray()

The code will look for PFAM domains in the database . If it can't find any it'll 
throw a value error.
Malformed partitions should be caught
"""

def test_partition_getpfamdomains_valid_pdb(mock_database):
    """
    Test that if a PFAM domain is found in the database, then
    the pfam domain array is correctly returned.
    """
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    array = partition.getPFAMDomainArray()
    for element in array[:3]:
        assert element == 2
    for element in array[3:]:
        assert element == 1
    


def test_partition_getpfamdomains_invalid_pdb(mock_database):
    """
    Test that if no PFAM domains are located in the database, then
    no array is returned and an exception is thrown.
    """
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '2vcr',
        'N': 10,
        'edgelistid': ObjectId('58dcf13fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    with pytest.raises(ValueError):
        array = partition.getPFAMDomainArray()
    
"""
tests for plotPymolStructure

We can't run PyMoL so we gotta test that the input script (the pml file) is correctly generated
and assume the rest is fine.

Inputs:
- a partition (atomic or residue), perhaps with multiple levels
- level
- outputpng

- What happens if the level is too high?
- What if we want a PNG?
- What if we dont?

"""


def mock_osremove(monkeypatch):
    """
    Patch os.remove so that the temporary PyMOL files are retained.
    
    I'll then need to clean these files after checking them manually.    
    """
    print("Remove functionality has been monkeypatched; temp files retained")
    # assert 0

def test_partition_plotpymolstructure_valid(mock_database, monkeypatch):
    """
    Test the default settings work successfully by checking the generating pml file.
    """
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    
    removeBackup = os.remove
    monkeypatch.setattr("os.remove", mock_osremove)
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    partition.plotPymolStructure()
    monkeypatch.setattr("os.remove", removeBackup)
    # Inspect the temp.pml file

    # Key points:
    # - as many levels as there are Infomap levels
    # - the correct number of residues
    # - the correct colouring
    
    expectedPml = """load temp.pdb, 1ubq_0
load temp.pdb, 1ubq_1
alter 1ubq_0 and (resi 18 or resi 19 or resi 20 or resi 21 or resi 22 or resi 23 or resi 24 or resi 51 or resi 52 or resi 53 or resi 54 or resi 55 or resi 56 or resi 57 or resi 58 or resi 59 or resi 60 or resi 61 ), b=0.16666666666666666
alter 1ubq_0 and (resi 36 or resi 37 or resi 38 or resi 39 or resi 40 or resi 41 or resi 42 or resi 69 or resi 70 or resi 71 or resi 72 or resi 73 or resi 74 or resi 75 or resi 76 ), b=0.3333333333333333
alter 1ubq_0 and (resi 1 or resi 2 or resi 3 or resi 4 or resi 14 or resi 15 or resi 16 or resi 17 or resi 62 or resi 63 or resi 64 or resi 65 ), b=0.5
alter 1ubq_0 and (resi 43 or resi 44 or resi 45 or resi 46 or resi 47 or resi 48 or resi 49 or resi 50 or resi 66 or resi 67 or resi 68 ), b=0.6666666666666666
alter 1ubq_0 and (resi 25 or resi 26 or resi 27 or resi 28 or resi 29 or resi 30 or resi 31 or resi 32 or resi 33 or resi 34 or resi 35 ), b=0.8333333333333334
alter 1ubq_0 and (resi 5 or resi 6 or resi 7 or resi 8 or resi 9 or resi 10 or resi 11 or resi 12 or resi 13 ), b=1.0
alter 1ubq_1 and (resi 59 ), b=0.013157894736842105
alter 1ubq_1 and (resi 23 ), b=0.02631578947368421
alter 1ubq_1 and (resi 56 ), b=0.039473684210526314
alter 1ubq_1 and (resi 61 ), b=0.05263157894736842
alter 1ubq_1 and (resi 21 ), b=0.06578947368421052
alter 1ubq_1 and (resi 55 ), b=0.07894736842105263
alter 1ubq_1 and (resi 54 ), b=0.09210526315789473
alter 1ubq_1 and (resi 22 ), b=0.10526315789473684
alter 1ubq_1 and (resi 24 ), b=0.11842105263157894
alter 1ubq_1 and (resi 19 ), b=0.13157894736842105
alter 1ubq_1 and (resi 52 ), b=0.14473684210526316
alter 1ubq_1 and (resi 18 ), b=0.15789473684210525
alter 1ubq_1 and (resi 58 ), b=0.17105263157894737
alter 1ubq_1 and (resi 51 ), b=0.18421052631578946
alter 1ubq_1 and (resi 57 ), b=0.19736842105263158
alter 1ubq_1 and (resi 60 ), b=0.21052631578947367
alter 1ubq_1 and (resi 20 ), b=0.2236842105263158
alter 1ubq_1 and (resi 53 ), b=0.23684210526315788
alter 1ubq_1 and (resi 41 ), b=0.25
alter 1ubq_1 and (resi 69 ), b=0.2631578947368421
alter 1ubq_1 and (resi 40 ), b=0.27631578947368424
alter 1ubq_1 and (resi 42 ), b=0.2894736842105263
alter 1ubq_1 and (resi 36 ), b=0.3026315789473684
alter 1ubq_1 and (resi 72 ), b=0.3157894736842105
alter 1ubq_1 and (resi 71 ), b=0.32894736842105265
alter 1ubq_1 and (resi 38 ), b=0.34210526315789475
alter 1ubq_1 and (resi 70 ), b=0.35526315789473684
alter 1ubq_1 and (resi 37 ), b=0.3684210526315789
alter 1ubq_1 and (resi 39 ), b=0.3815789473684211
alter 1ubq_1 and (resi 73 ), b=0.39473684210526316
alter 1ubq_1 and (resi 74 ), b=0.40789473684210525
alter 1ubq_1 and (resi 75 ), b=0.42105263157894735
alter 1ubq_1 and (resi 76 ), b=0.4342105263157895
alter 1ubq_1 and (resi 4 ), b=0.4473684210526316
alter 1ubq_1 and (resi 3 ), b=0.4605263157894737
alter 1ubq_1 and (resi 15 ), b=0.47368421052631576
alter 1ubq_1 and (resi 1 ), b=0.4868421052631579
alter 1ubq_1 and (resi 2 ), b=0.5
alter 1ubq_1 and (resi 17 ), b=0.5131578947368421
alter 1ubq_1 and (resi 65 ), b=0.5263157894736842
alter 1ubq_1 and (resi 64 ), b=0.5394736842105263
alter 1ubq_1 and (resi 14 ), b=0.5526315789473685
alter 1ubq_1 and (resi 63 ), b=0.5657894736842105
alter 1ubq_1 and (resi 62 ), b=0.5789473684210527
alter 1ubq_1 and (resi 16 ), b=0.5921052631578947
alter 1ubq_1 and (resi 45 ), b=0.6052631578947368
alter 1ubq_1 and (resi 43 ), b=0.618421052631579
alter 1ubq_1 and (resi 67 ), b=0.631578947368421
alter 1ubq_1 and (resi 44 ), b=0.6447368421052632
alter 1ubq_1 and (resi 50 ), b=0.6578947368421053
alter 1ubq_1 and (resi 68 ), b=0.6710526315789473
alter 1ubq_1 and (resi 66 ), b=0.6842105263157895
alter 1ubq_1 and (resi 49 ), b=0.6973684210526315
alter 1ubq_1 and (resi 48 ), b=0.7105263157894737
alter 1ubq_1 and (resi 46 ), b=0.7236842105263158
alter 1ubq_1 and (resi 47 ), b=0.7368421052631579
alter 1ubq_1 and (resi 27 ), b=0.75
alter 1ubq_1 and (resi 30 ), b=0.7631578947368421
alter 1ubq_1 and (resi 31 ), b=0.7763157894736842
alter 1ubq_1 and (resi 26 ), b=0.7894736842105263
alter 1ubq_1 and (resi 29 ), b=0.8026315789473685
alter 1ubq_1 and (resi 25 ), b=0.8157894736842105
alter 1ubq_1 and (resi 33 ), b=0.8289473684210527
alter 1ubq_1 and (resi 34 ), b=0.8421052631578947
alter 1ubq_1 and (resi 28 ), b=0.8552631578947368
alter 1ubq_1 and (resi 32 ), b=0.868421052631579
alter 1ubq_1 and (resi 35 ), b=0.881578947368421
alter 1ubq_1 and (resi 5 ), b=0.8947368421052632
alter 1ubq_1 and (resi 13 ), b=0.9078947368421053
alter 1ubq_1 and (resi 6 ), b=0.9210526315789473
alter 1ubq_1 and (resi 7 ), b=0.9342105263157895
alter 1ubq_1 and (resi 12 ), b=0.9473684210526315
alter 1ubq_1 and (resi 11 ), b=0.9605263157894737
alter 1ubq_1 and (resi 8 ), b=0.9736842105263158
alter 1ubq_1 and (resi 9 ), b=0.9868421052631579
alter 1ubq_1 and (resi 10 ), b=1.0
#formatting
bg_color white
hide all
#show sticks
show cartoon
spectrum b, rainbow,  minimum=0, maximum=1
set opaque_background=0
set antialias = on
set line_smooth = 1
set depth_cue = 1
set specular = 1
set surface_quality = 1
set stick_quality = 15
set sphere_quality = 2
set ray_trace_fog = 0.8
set light = (-0.2,0,-1)

set ray_shadows, 0
set surface_mode, 1
set cartoon_side_chain_helper,on
rebuild
save 1ubq.pse 
"""   
    with open("temp.pml") as flines:
        generatedPml = flines.read()

    assert generatedPml == expectedPml
    
    # Tidy
    os.remove("temp.pdb")
    os.remove("temp.pml")
    # assert 0

def test_partition_plotpymolstructure_valid_level(mock_database, monkeypatch):
    """
    Test the default settings work successfully by checking the generating pml file.
    """
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    
    removeBackup = os.remove
    monkeypatch.setattr("os.remove", mock_osremove)
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    partition.plotPymolStructure(level=0)
    monkeypatch.setattr("os.remove", removeBackup)
    # Inspect the temp.pml file

    # Key points:
    # - as many levels as there are Infomap levels
    # - the correct number of residues
    # - the correct colouring
    
    expectedPml = """load temp.pdb, 1ubq_0
alter 1ubq_0 and (resi 18 or resi 19 or resi 20 or resi 21 or resi 22 or resi 23 or resi 24 or resi 51 or resi 52 or resi 53 or resi 54 or resi 55 or resi 56 or resi 57 or resi 58 or resi 59 or resi 60 or resi 61 ), b=0.16666666666666666
alter 1ubq_0 and (resi 36 or resi 37 or resi 38 or resi 39 or resi 40 or resi 41 or resi 42 or resi 69 or resi 70 or resi 71 or resi 72 or resi 73 or resi 74 or resi 75 or resi 76 ), b=0.3333333333333333
alter 1ubq_0 and (resi 1 or resi 2 or resi 3 or resi 4 or resi 14 or resi 15 or resi 16 or resi 17 or resi 62 or resi 63 or resi 64 or resi 65 ), b=0.5
alter 1ubq_0 and (resi 43 or resi 44 or resi 45 or resi 46 or resi 47 or resi 48 or resi 49 or resi 50 or resi 66 or resi 67 or resi 68 ), b=0.6666666666666666
alter 1ubq_0 and (resi 25 or resi 26 or resi 27 or resi 28 or resi 29 or resi 30 or resi 31 or resi 32 or resi 33 or resi 34 or resi 35 ), b=0.8333333333333334
alter 1ubq_0 and (resi 5 or resi 6 or resi 7 or resi 8 or resi 9 or resi 10 or resi 11 or resi 12 or resi 13 ), b=1.0
#formatting
bg_color white
hide all
#show sticks
show cartoon
spectrum b, rainbow,  minimum=0, maximum=1
set opaque_background=0
set antialias = on
set line_smooth = 1
set depth_cue = 1
set specular = 1
set surface_quality = 1
set stick_quality = 15
set sphere_quality = 2
set ray_trace_fog = 0.8
set light = (-0.2,0,-1)

set ray_shadows, 0
set surface_mode, 1
set cartoon_side_chain_helper,on
rebuild
save 1ubq.pse 
"""   
    with open("temp.pml") as flines:
        generatedPml = flines.read()

    assert generatedPml == expectedPml
    
    # Tidy
    os.remove("temp.pdb")
    os.remove("temp.pml")
    # assert 0


def test_partition_plotpymolstructure_invalid_level(mock_database, monkeypatch):
    """
    Test that invalid level choices (i.e. ints <0 or >1 in this case, or any non-int)
    are rejected correctly.
    """
    db = proteinnetworks.database.Database(password="bla")
    partitionArgs = {
        'pdbref': '1ubq',
        'N': 10,
        'edgelistid': ObjectId('58dbe03fef677d54224a01da'),
        'detectionmethod': 'Infomap',
        'r': -1,
        "database": db
    }
    
    removeBackup = os.remove
    monkeypatch.setattr("os.remove", mock_osremove)
    partition = proteinnetworks.partition.Partition(**partitionArgs)
    with pytest.raises(IndexError):
        partition.plotPymolStructure(level=-10)

    # assert 0

    