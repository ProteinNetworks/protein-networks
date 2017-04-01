"""
Create sets of indices which enforce uniqueness for the given parameter set.

I don't want multiple copies of an document in my database (for space and
also because it implies something has gone terribly wrong.) We can fix this
by creating a compound index which is unique, such that a certain combination of
parameters cannot be repeated.

But different documents have different fields.
So we want partial, compound, unique indices covering all documents.

Four unique partial indices:
(pdbref) if doctype==pdbfile
(pdbref, edgelisttype, hydrogenstatus, scaling) if doctype == edgelist
(pdbref. edgelistid, detectionmethod, r) if detectionmethod == "AFG"
(pdbref, edgelistid, detectionmethod) if detectionmethod == "Infomap
"""
import pymongo

if __name__ == "__main__":
    password = input("password: ").strip()
    client = pymongo.MongoClient(
        "mongodb://writeAccess:" + password + "@127.0.0.1/proteinnetworks",
        serverSelectionTimeoutMS=1000)
    db = client.proteinnetworks
    collection = db.proteinnetworks

    collection.drop_indexes()
    # print(collection.find({"doctype": "pdbfile"}).count())
    collection.create_index(
        "pdbref", unique=True, partialFilterExpression={"doctype": "pdbfile"})

    collection.create_index(
        [("pdbref", pymongo.ASCENDING), ("edgelisttype", pymongo.ASCENDING),
         ("hydrogenstatus", pymongo.ASCENDING),
         ("scaling", pymongo.ASCENDING)],
        unique=True,
        partialFilterExpression={"doctype": "edgelist"})
    collection.create_index(
        [("pdbref", pymongo.ASCENDING), ("edgelistid", pymongo.ASCENDING),
         ("detectionmethod", pymongo.ASCENDING), ("r", pymongo.ASCENDING)],
        unique=True,
        partialFilterExpression={"detectionmethod": "AFG"})
    collection.create_index(
        [("pdbref", pymongo.ASCENDING), ("edgelistid", pymongo.ASCENDING),
         ("detectionmethod", pymongo.ASCENDING), ("N", pymongo.ASCENDING)],
        unique=True,
        partialFilterExpression={"detectionmethod": "Infomap"})
