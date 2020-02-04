import numpy as np
import graph_tool.all as gt
from collections import defaultdict

def buildGraph(edgefilename, nodenamefile=None):
    g = gt.Graph(directed=True)
    regulation = g.new_edge_property('int')

    names = []
    if nodenamefile!=None:
        with open(nodenamefile, 'r') as namelist:
            for line in namelist:
                data = line.split()
                names.append(data[0])

    if nodenamefile==None:
        with open(edgefilename, 'r') as edgelist:
            for line in edgelist:
                data = line.split()

                cur_edge = g.add_edge(int(data[0]), int(data[1]), add_missing=True)
                if data[2]=='positive':
                    regulation[cur_edge] = 0
                elif data[2]=='negative':
                    regulation[cur_edge] = 1
                elif data[2]=='dual':
                    regulation[cur_edge] = 2
    else:
        nodename = g.new_vertex_property('string')
        with open(edgefilename, 'r') as edgelist:
            for line in edgelist:
                data = line.split()

                cur_edge = g.add_edge(int(data[0]), int(data[1]), add_missing=True)
                nodename[int(data[0])] = names[int(data[0])]
                nodename[int(data[1])] = names[int(data[1])]
                if data[2]=='positive':
                    regulation[cur_edge] = 0
                elif data[2]=='negative':
                    regulation[cur_edge] = 1
                elif data[2]=='dual':
                    regulation[cur_edge] = 2
        g.vertex_properties['node_names'] = nodename

    g.edge_properties['regulation'] = regulation
    return g


def edgefromSet_optimal(arr, graph, pivot, pivotnode_to_index, regulation):
    ''' given a pivot set and an 'number of received 
        information' matrix 'arr', with size (n_edgetype, 
        len(pivot_sucessor)), calculates the value of the entries
        of 'arr'.   '''

    pivot_nodes = pivot.get_nodes()
    for node in pivot_nodes:
        out_edges = graph.get_out_edges(node, [graph.edge_index])
        for out in out_edges:
            reg = regulation[out[2]]
            correct_index = pivotnode_to_index[out[1]]
            arr[reg][correct_index] += 1


# Input-tree stability
def edgefromSet(arr, graph, refinement_set, regulation):
    set_nodes = refinement_set.get_nodes()
    for setnode in set_nodes:
        out_edges = graph.get_out_edges(setnode, [graph.edge_index])
        for edge in out_edges:
            reg = regulation[edge[2]]
            arr[reg][edge[1]] += 1

# Output-tree stability
def edgetoSet(arr, graph, refinement_set, regulation):
    set_nodes = refinement_set.get_nodes()
    for setnode in set_nodes:
        in_edges = graph.get_in_edges(setnode, [graph.edge_index])
        for edge in in_edges:
            reg = regulation[edge[2]]
            arr[reg][edge[0]] += 1

# Input solitaire
def IDENTIFY_SOLITAIRE(graph, node):
    in_neighbors = graph.get_in_neighbors(node)
    n_in = in_neighbors.shape[0]
    if n_in==0: return 0    # Solitaire and do not feeds itself.
    #elif node in in_neighbors: return 1
    elif np.all(in_neighbors==node): return 1 # Solitaire but feeds itself.
    else: return -1 # Not solitaire.

# Output solitaire
def OUT_IDENTIFY_SOLITAIRE(graph, node):
    out_neighbors = graph.get_out_neighbors(node)
    n_out = out_neighbors.shape[0]
    if n_out==0: return 0    # Solitaire and do not feeds itself.
    elif np.all(out_neighbors==node): return 1 # Solitaire but feeds itself.
    else: return -1 # Not solitaire.

######################################################################
class VisitedinBFS(gt.BFSVisitor):
    def __init__(self):
        self.visitedlist = []

    def discover_vertex(self, u):
        self.visitedlist.append(int(u))
        
def GetStrongCCompOfNode(graph, vertex):
    # Extract nodes from a bfs search starting at 'vertex'.
    nodes = VisitedinBFS()
    gt.bfs_search(graph, source=graph.vertex(vertex), visitor=nodes)

    # For each one of the extracted nodes, we perform a bfs search 
    # starting at this node. If 'vertex' is reached, then both belong
    # to the same strongly connected component.
    stronglycc = []
    for w in nodes.visitedlist:
        aux = VisitedinBFS()
        gt.bfs_search(graph, source=graph.vertex(w), visitor=aux)
        if vertex in aux.visitedlist:
            stronglycc.append(w)
    return stronglycc        

###################### GET INFORMATION UTILS #########################
def GetNumberFibers(partition):
    count = 0
    nodes_in_fiber = 0
    for block in partition:
        if block.get_number_nodes()>1:
            count +=1 
            nodes_in_fiber += block.get_number_nodes()
    return (count, nodes_in_fiber)

def PrintFibers(partition, graph, name=False):
    if name:
        for block in partition:
            block.show_nodes_name(graph)
    else:
        for block in partition:
            block.show_nodes()
#######################################################################
