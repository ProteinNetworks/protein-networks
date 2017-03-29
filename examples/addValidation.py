"""
Add document validation rules to the MongoDB collection.

Validation rules:

    "pdbref" exists AND
    "data" exists AND
        "doctype" is "edgelist" AND
            ("edgelisttype" is "atomic" OR "residue") AND
            ("hydrogenstatus" is "noH" OR "Hatoms" OR "Hbonds") AND
            ("scaling" is a float)
        OR
        "doctype" is "partition" AND
            ("edgelistid" is an ObjectID) AND
            (("detectionmethod" is "AFG" AND "r" is a float) OR
             ("detectionmethod" is "Infomap")
       OR
       "doctype" is "pdbfile"


NB this might be too large for MongoDB to handle (>16MB)
"""
import pymongo
import json

if __name__ == "__main__":
    client = pymongo.MongoClient()
    db = client.proteinnetworks
    collection = db.proteinnetworks
    validator = {
        "$and": [{
            "pdbref": {
                "$type": "string"
            }
        }, {
            "data": {
                "$exists": True
            }
        }, {
            "$or": [{
                "$and": [{
                    "doctype": "edgelist"
                }, {
                    "$or": [{
                        "edgelisttype": "atomic"
                    }, {
                        "edgelisttype": "residue"
                    }]
                }, {
                    "$or": [{
                        "hydrogenstatus": "noH"
                    }, {
                        "hydrogenstatus": "Hatoms"
                    }, {
                        "hydrogenstatus": "Hbonds"
                    }]
                }, {
                    "scaling": {
                        "$type": "number"
                    }
                }]
            }, {
                "$and": [{
                    "doctype": "partition"
                }, {
                    "edgelistid": {
                        "$type": "objectId"
                    }
                }, {
                    "$or": [{
                        "$and": [{
                            "detectionmethod": "AFG"
                        }, {
                            "r": {
                                "$type": "number"
                            }
                        }]
                    }, {
                        "detectionmethod": "Infomap"
                    }]
                }]
            }, {
                "doctype": "pdbfile"
            }]
        }]
    }
    db.command("collMod", "proteinnetworks", validator=validator)
    collections = db.command("listCollections")
    print(json.dumps(collections, indent=2))
