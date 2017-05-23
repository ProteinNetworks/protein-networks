"""
Adds the pdb-pfam mappings to the database.

Example lines in the mapping file:
PDB_ID	CHAIN_ID	PdbResNumStart	PdbResNumEnd	PFAM_ACC	PFAM_Name	PFAM_desc	eValue
1AHU	A	72	213	PF01565.19	FAD_binding_4	FAD binding domain	8.2E-26

Import as:
{
    "pdbref": bla,
    "doctype": mapping,
    "mappingtype": PFAM,
    "data": { "chainid": A, "startresidue": 1, "endresidue": 10, "pfamref": bla}
}
"""
import pymongo
if __name__ == "__main__":
    password = input("password: ").strip()
    client = pymongo.MongoClient(
        "mongodb://writeAccess:" + password + "@127.0.0.1/proteinnetworks",
        serverSelectionTimeoutMS=1000)
    db = client.proteinnetworks
    collection = db.proteinnetworks

    with open("/home/will/MainProject/pdb_pfam_mapping.txt") as flines:
        for line in flines:
            lineContents = line.strip().split("\t")
            pdbRef = lineContents[0]
            chainId = lineContents[1]
            startRes = lineContents[2]
            endRes = lineContents[3]
            pfamRef = lineContents[4]
            document = {
                "pdbref": pdbRef.lower(),
                "doctype": "mapping",
                "mappingtype": "PFAM",
                "data": {
                    "chainid": chainId,
                    "startresidue": startRes,
                    "endresidue": endRes,
                    "pfamref": pfamRef
                }
            }
            result = collection.insert_one(document)
            print(result.inserted_id)
