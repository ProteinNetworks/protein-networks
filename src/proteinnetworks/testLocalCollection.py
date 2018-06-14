import proteinnetworks
from bson import ObjectId
import json
from pprint import pprint
db = proteinnetworks.database.Database(local=True)

print("Count:",db.collection.count())


depositionArgs = {
    'pdbref': '2vc5',
    'edges':
    [[2, 1, 44], [3, 1, 40], [3, 2, 56], [4, 2, 56], [4, 3, 70], [5, 3, 23]],
    'edgelisttype': 'residue',
    'hydrogenstatus': 'noH',
    'scaling': 4.5
}
edgelistId = db.depositEdgelist(**depositionArgs)
print("Count:",db.collection.count())
depositionArgs = {
    'pdbref': '2vc5',
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
    'edgelistid': ObjectId(edgelistId),
    'detectionmethod': 'Infomap',
    'r': -1
}
resultId = db.depositPartition(**depositionArgs)
print("Count:",db.collection.count())

query = {
    'pdbref': '2vc5',
    'N': 10,
    'edgelistid': ObjectId(edgelistId),
    'detectionmethod': 'Infomap',
}
cursor = db.collection.find(query)

print(cursor)
db.depositPartition(**depositionArgs)
