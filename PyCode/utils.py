import numpy as np
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


# Much more efficient
def edgefromSet(arr, graph, refinement_set, regulation):
    set_nodes = refinement_set.get_nodes()
    for setnode in set_nodes:
        out_edges = graph.get_out_edges(setnode, [graph.edge_index])
        for edge in out_edges:
            reg = regulation[edge[2]]
            arr[reg][edge[1]] += 1


def IDENTIFY_SOLITAIRE(graph, node):
    in_neighbors = graph.get_in_neighbors(node)
    n_in = in_neighbors.shape[0]
    if n_in==0: return 0    # Solitaire and do not feeds itself.
    elif node in in_neighbors and n_in==1: return 1 # Solitaire but feeds itself.
    else: return -1 # Not solitaire.

def GetNumberFibers(partition):
    count = 0
    nodes_in_fiber = 0
    for block in partition:
        if block.get_number_nodes()>1:
            count +=1 
            nodes_in_fiber += block.get_number_nodes()
    return (count, nodes_in_fiber)

def PrintFibers(partition):
    for block in partition:
        block.show_nodes()
