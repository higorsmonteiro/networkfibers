import graph_tool.all as gt

def buildGraph(edgefilename):
    g = gt.Graph(directed=True)
    regulation = g.new_edge_property('int')

    with open(edgefilename) as edgelist:
        for line in edgelist:
            data = line.split()
            
            cur_edge = g.add_edge(int(data[0]), int(data[1]), add_missing=True)
            if data[2]=='positive':
                regulation[cur_edge] = 0
            elif data[2]=='negative':
                regulation[cur_edge] = 1
            elif data[2]=='dual':
                regulation[cur_edge] = 2
    g.edge_properties['regulation'] = regulation
    return g


def edgefromSet(arr, graph, node, refinement_set, regulation, reg_type):
    '''
        Given a refinement set and an edge type, the function calculates
        the number of edges of that type that 'node' receives from the 
        refinement set, and stores it in the 'arr'.
    '''
    in_edges = graph.get_in_edges(node, [graph.edge_index])
    ref_set = refinement_set.get_nodes()

    for edge in in_edges:
        if edge[0] in ref_set and regulation[edge[2]]==reg_type:
            arr[node] += 1


def IDENTIFY_SOLITAIRE(graph, node):
    in_neighbors = graph.get_in_neighbors(node)
    n_in = in_neighbors.shape[0]
    if n_in==0: return 0    # Solitaire and do not feeds itself.
    elif node in in_neighbors and n_in==1: return 1 # Solitaire but feeds itself.
    else: return -1 # Not solitaire.
