import graph_tool.all as gt

def buildGraph(edgefilename):
    g = gt.Graph(directed=True)
    regulation = []

    with open(edgefilename) as edgelist:
        for line in edgelist:
            data = line.split("\t")
            regulation.append(data[2])
            g.add_edge(int(data[0]), int(data[1]), add_missing=True)
    return g

def edgefromSet(arr_fromSet, graph, node, refinement_set, reg_type):
    in_neighbors = graph.get_in_neighbors(node)

    # and regulation type?
    for neigh in in_neighbors:
        if neigh in refinement_set:
            arr_fromSet[node] += 1


def IDENTIFY_SOLITAIRE(graph, node):
    in_neighbors = graph.get_in_neighbors(node)
    n_in = in_neighbors.shape[0]
    if n_in==0: return 0    # Solitaire and do not feeds itself.
    elif node in in_neighbors and n_in==1: return 1 # Solitaire but feeds itself.
    else: return -1 # Not solitaire.
