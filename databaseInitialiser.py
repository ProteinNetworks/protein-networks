"""
Creates the MongoDB database, from all the data I have scattered around my filesystem.

MongoDB database format:
database: proteinnetworks
collection: proteinnetworks
document has fields:
    _id: bla
    pdbref: bla
    doctype: (edgelist | partition | pdbfile ) currently

if doctype == edgelist:
    edgelisttype: (atomic | residue)
    hydrogenstatus: (noH | Hatoms | Hbonds) (i.e. preprocessed with nothing, HAAD, or Stride)
    scaling: the value of the cutoff used
    data: the actual edgelist

if doctype == partition:
    edgelistused: The _id of the edgelist used in generating the partition
    detectionmethod: (AFG | Infomap) as it stands

    if detectionmethod == AFG:
        r : The AFG parameter used in generating the partition
        data: The 1D array giving the partition

    if detectionmethod == Infomap:
        data: A 2D numpy array giving the data

if doctype == pdbfile:
    data: The PDBfile itself


NB this might be too large for MongoDB to handle (>16MB)
"""


if __name__ == "__main__":
    print("Hello World")
