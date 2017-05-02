# Run the SBM on all edgelists
#
import graph_tool
import graph_tool.inference
import os
import proteinnetworks

db = proteinnetworks.database.Database()

edgelists = db.collection.find({
    "doctype": "edgelist",
    "edgelisttype": "residue",
    "scaling": 4.0
})

for edgelist in edgelists:
    edgelistFilename = "Data/{}.{}.dat".format(edgelist['pdbref'],
                                               edgelist['scaling'])
    if os.path.exists(edgelistFilename):
        continue
    with open(edgelistFilename, mode='w') as flines:
        flines.write(
            "\n".join(" ".join(map(str, x)) for x in edgelist['data']))

    g = graph_tool.load_graph_from_csv(
        edgelistFilename,
        directed=False,
        skip_first=False,
        csv_options={"delimiter": " "})

    state = graph_tool.inference.minimize_nested_blockmodel_dl(
        g, deg_corr=True)
    hierarchical_partition = []
    for i, innerstate in enumerate(state.get_levels()):
        hierarchical_partition.append(
            list(state.project_level(i).b.get_array()))
    hierarchical_partition = hierarchical_partition[::-1][1:]

    SBMFilename = "Data/{}.{}.sbm".format(edgelist['pdbref'],
                                          edgelist['scaling'])
    with open(SBMFilename, mode='w') as flines:
        flines.write(
            "\n".join(" ".join(map(str, x)) for x in hierarchical_partition))
