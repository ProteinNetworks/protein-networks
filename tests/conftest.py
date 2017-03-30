"""Mock for the MongoDB database, used in all tests."""
import pytest


@pytest.fixture(autouse=True)
def mock_database(monkeypatch):
    """Monkeypatch the pymongo.client() and adds some test data."""
    class Gary:
        """The mock MongoClient. Does nothing, simply returns."""

        @classmethod
        def __init__(self, inputString, serverSelectionTimeoutMS):
            pass

        class admin:
            """
            The admin.command which is used solely to test the connection.

            Again, just return (for now)
            """

            def command(inputString):
                print("Gary used instead of database")

            pass

        class proteinnetworks:
            """The proteinnetworks "database"."""

            class proteinnetworks:
                """
                The proteinnetworks "collection".

                This is where, in the real system, the data is stored.
                """

                def find(query):
                    """Pretend to extract a protein."""
                    class Cursor:

                        def __init__(self):
                            self.doc = {"data": 5, "_id": 10}
                            pass

                        def count(self):
                            return 1

                        def __getitem__(self, i):
                            return self.doc

                    return Cursor()

    monkeypatch.setattr("pymongo.MongoClient", Gary)
