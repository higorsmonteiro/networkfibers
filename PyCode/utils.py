import numpy as np
import graph_tool.all as gt
from collections import defaultdict

def defineGraph(edgefilename, nodenamefile=None):
    g = gt.Graph(directed=True)
    fiber_index = g.new_vertex_property('int')
    information_type = g.new_edge_property('int')
    if nodenamefile!=None:
        nodename =  g.new_vertex_property('string')

    # Counts the # of different types of information in the file.
    edgetype_set = set()
    with open(edgefilename, 'r') as edgelist:
        for line in edgelist:
            data = line.split()
            # The 3-th column must be the edge type.
            edgetype_set.add(data[2])
    
    # For each type string associate it to an integer index.
    count = 0
    edgetype_to_index = defaultdict(lambda:-1)
    for info_type in edgetype_set:
        edgetype_to_index[info_type] = count
        count+=1
    ##############################################################
    
    ''' If a file name containing strings to each node
        is provided, then we store them in 'names' and
        then we define the network. Otherwise, we build
        the network only using the edgelist file.  '''
    names = []
    if nodenamefile!=None:
        with open(nodenamefile, 'r') as namelist:
            for line in namelist:
                names.append(line.split()[0])

        with open(edgefilename, 'r') as edgelist:
            for line in edgelist:
                data = line.split()
                sender = int(data[0])
                receiver = int(data[1])

                # 'add_missing' True because sender and receiver may not exists in the network yet.
                edge = g.add_edge(sender, receiver, add_missing=True)
                nodename[sender] = names[sender]
                nodename[receiver] = names[receiver]
                # use the dict 'edgetype_to_index' to give the its integer index.
                information_type[edge] = edgetype_to_index[data[2]]
        g.vertex_properties['nodename'] = nodename
        ###################################################
    else:
        with open(edgefilename, 'r') as edgelist:
            for line in edgelist:
                data = line.split()
                sender = int(data[0])
                receiver = int(data[1])

                edge = g.add_edge(sender, receiver, add_missing=True)
                information_type[edge] = edgetype_to_index[data[2]]
    
    for v in g.get_vertices(): fiber_index[v] = -1
    g.vertex_properties['fiber_index'] = fiber_index
    g.edge_properties['regulation'] = information_type
    return g, len(edgetype_set)


def calc_R(R, graph, pivot, f, regulation):
    ''' given a pivot set and an 'number of received 
        information' matrix 'R', with size (n_edgetype, 
        len(pivot_sucessor)), calculates the value of the entries
        of 'R'.   '''

    pivot_nodes = pivot.get_nodes()
    for node in pivot_nodes:
        out_edges = graph.get_out_edges(node, [graph.edge_index])
        for out in out_edges:
            reg = regulation[out[2]]
            correct_index = f[out[1]]
            R[reg,correct_index] += 1

def is_unstable(arr_2d):
    ''' For a matrix to be stable, for each row all the
        columns must be equal. Otherwise, the matrix is
        unstable and this function return True.   '''
    for row in arr_2d:
        if not np.all(row==row[0]):
            return True # The matrix is unstable.
    return False

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
