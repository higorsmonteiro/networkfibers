import numpy as np
from fiber import *
import graph_tool.all as gt
from collections import Counter, defaultdict

def Initialization(graph):
    ''' 
        Separate each node in its corresponding SCC using
        the 'label_components' function from graph_tool
        library.
    '''
    N = graph.get_vertices().shape[0]
    label_scc, hist = gt.label_components(graph, directed=True)

    fibers_listing = defaultdict(list)
    for v in graph.get_vertices():
        label = label_scc[v]
        fibers_listing[label].append(int(v))
    
    ''' 
        'fibers' now contains for each SCC label, 
        all nodes belonging to it.
    '''

    scc = []
    N_scc = hist.shape[0]
    # Insert each node in its correct SCC object.
    for scc_j in range(N_scc):
        scc.append(StrongComponent())
        nodes_list = fibers_listing[scc_j]
        for n in nodes_list:
            scc[scc_j].insert_node(n)
    '''
        'scc[j]' is an object containing the information
        about the nodes inside the j-th SCC.
    '''

    ''' Defines if each SCC receives or not input
        from other components not itself. '''
    for strong in scc:
        strong.check_input(graph)
    
    # Count the number of colors.
    color_index = 0
    fibers = [FiberBlock()]
    fibers[0].index = color_index
    for strong in scc:
        input_bool = strong.get_input_bool()
        if input_bool==True:
            fibers[0].insert_nodelist(strong.nodes)
        else:
            color_index += 1
            fibers.append(FiberBlock())
            fibers[-1].index = color_index
            fibers[-1].insert_nodelist(strong.nodes)

    return fibers

def set_colors(graph, fibers):
    node_colors = graph.vp.node_colors
    
    for index, fiberblock in enumerate(fibers):
        for v in fiberblock.fibernodes:
            node_colors[v] = index

# Define the ISCV for each node. At this step, each ISCV has size 'ncolor'.
def set_ISCV(graph, ncolor):
    # intrinsic properties
    iscv = graph.vp.iscv
    node_colors = graph.vp.node_colors

    # initializate iscv
    for v in graph.get_vertices():
        iscv[v] = ""

    # define the iscv for each node.
    for v in graph.get_vertices():
        in_neighbors = graph.get_in_neighbors(v)
        input_colors = []
        for neigh in in_neighbors:
            input_colors.append(node_colors[neigh])
        # Counts how many inputs the 'node' receives from each color.
        Colors_counter = Counter(input_colors)
        for k in range(ncolor):
            # 'Color_counter[k] returns the number of inputs from color 'k'.
            iscv[v] += str(Colors_counter[k])

def split_fiber(fibernodes, fiber_iscv):
    


def print_colors(graph):
    node_colors = graph.vp.node_colors
    for v in graph.get_vertices():
        print(v, node_colors[v])

def number_colors(graph):
    node_colors = list(graph.vp.node_colors.a)
    print(len(set(node_colors)))
