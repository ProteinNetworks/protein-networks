from bson import ObjectId


class LocalCollection:
    """
    A very limited in-memory "database" to be used if no MongoDB instance can be found.

    Stores records as a list of dicts, manipulated with the following methods:
    - collection.find():
        given a set of parameters (including wildcards such as $exists etc) as a dict,
        return a list of the dicts in the database matching this description.

    - collection.find_one():
        as above, but only return one. Used in this code to pull things by their ObjectId.

    - collection.insert_one():
        given a dict, add an ObjectId, push the record, return the id.

    - count():
        return the number of records in the db.
    """

    def __init__(self):
        """Initialise the empty list of dicts."""
        self.storageList = [{"pdbref": "1ubq", "data": "bla"}, {"pdbref": "2ubq", "data": "bla2"}]

    def find(self, query):
        """
        Return a 'cursor' which behaves like a generator with a count method.

        """
        class Cursor(list):
            """Extend the list class with a count method that does the same thing as len()."""
            def count(self):
                return len(self)

        subset = []
        for record in self.storageList:
            for key, value in query.items():
                if type(value) == dict and "$exists" in value:
                    exists = value["$exists"]
                    # match if "exists" is False and key isn't in the record
                    # or if "exists" is True and key is in the record
                    match = (exists and key in record) or (not exists and (key not in record))
                    if not match:
                        break
                else:
                    if record[key] != value:
                        break
            else:
                subset.append(record)

        results = Cursor(subset)
        return results

    def find_one(self, query):
        """
        Return a single record.

        As we expect the local DB to be small, this can just be find() with a length check.
        """
        results = self.find(query)
        assert len(results) < 2
        return results

    def insert_one(self, record):
        """
        Push a dictionary to the "database", adding a BSON ObjectId, and return a Result
        (with an inserted_id attribute).
        """
        class Result:
            """A container for the inserted_id, necessary to match the pymongo collection."""
            def __init__(self, id):
                self.inserted_id = id

        record["_id"] = ObjectId()
        self.storageList.append(record)
        result = Result(record["_id"])
        return result

    def count(self):
        return len(self.storageList)


collection = LocalCollection()
query = {"pdbref": "2ubq"}
results = collection.find(query)
results = collection.find_one(query)
print(results)
