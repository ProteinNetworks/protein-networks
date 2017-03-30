"""Mock for the MongoDB database, used in all tests."""
import pytest


@pytest.fixture(autouse=True)
def mock_database(monkeypatch):
    """Monkeypatch the pymongo.client() and adds some test data."""
    class Gary:
        @classmethod
        def __init__(self):
            print("Gary used instead of database")
            pass
    monkeypatch.setattr("proteinnetworks.database.Database", Gary)
