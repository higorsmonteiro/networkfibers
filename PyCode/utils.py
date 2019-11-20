import graph_tool.all as gt

def buildGraph(edgefilename):
    g = gt.Graph(directed=True)
    regulation = g.new_edge_property('int')

    with open(edgefilename) as edgelist:
        for line in edgelist:
            data = line.split()
            
            cur_edge = g.add_edge(int(data[0]), int(data[1]), add_missing=True)
            if data[2]=='positive\n':
                regulation[cur_edge] = 0
            elif data[2]=='negative\n':
                regulation[cur_edge] = 1
            elif data[2]=='dual\n':
                regulation[cur_edge] = 2
    return (g, regulation)

def edgefromSet(arr_fromSet, graph, node, refinement_set, reg_type):
    in_neighbors = graph.get_in_neighbors(node)

    # and regulation type?
    for neigh in in_neighbors:
        if neigh in refinement_set.get_nodes():
            arr_fromSet[node] += 1


def IDENTIFY_SOLITAIRE(graph, node):
    in_neighbors = graph.get_in_neighbors(node)
    n_in = in_neighbors.shape[0]
    if n_in==0: return 0    # Solitaire and do not feeds itself.
    elif node in in_neighbors and n_in==1: return 1 # Solitaire but feeds itself.
    else: return -1 # Not solitaire.
