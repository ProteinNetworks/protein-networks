# protein-networks
This repo stores the source code, any required database initialisation/configuration scripts, and the associated tests, for the networks section of my PhD (which can be found [here](https://wpg.io/phdThesis.pdf)).

[![Build Status](https://travis-ci.com/wllgrnt/protein-networks.svg?token=SQysmopgMsk5Ksy6iu6f&branch=master)](https://travis-ci.com/wllgrnt/protein-networks)

This code assumes that:

- There is a MongoDB instance in the default place. This is currently hardcoded in src/database.py. 
- The "Infomap" community detection method can be found in $PATH. This code predates the Python API, so needs the compiled version [here](https://github.com/mapequation/infomap). Code has been tested with version 0.18.5.

Todo:
- [ ] Add MongoDB configuration (even if hardcoded somewhere)
- [ ] Check test coverage
- [ ] Add sphinx docs 
- [ ] Add example notebooks

## Installation

Run `pip install . -I` from the repo root directory.

## Testing

Run `pytest` from the repo root directory.

## Usage

The code has three main parts:

`database.py`: Stores the Database class, which provides a wrapper to the MongoDB database, and handles fetching, storing, and validating data (see "scripts" for the validation schema used). The data can be stored locally using `local=True` on instance creation.

`network.py`: Stores the Network class, which handles network generation, plotting, etc. A network is generated with a given set of parameters, as follows:
```
    db = proteinnetworks.database.Database(password="bla")
    inputArgs = {
        "scaling": 4.5,
        "edgelisttype": "residue",
        "hydrogenstatus": "noH",
        "pdbref": "2vcr",
        "database": db
    }
    pn = proteinnetworks.network.Network(**inputArgs)
```

The database is queried for the matching network - if one is found in the database, it's returned. If not, the network is generated and stored in the database.

`partition.py`: Stores the Partition class, which uses Infomap to generate the community structure for a given network. As in the Network case, the database is first queried for a matching community structure, and if not found then Infomap is run, and the results stored.

`insight.py` stores convenience functions and classes for scoring partitions etc.

