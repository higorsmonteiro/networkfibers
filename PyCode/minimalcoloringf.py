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
    fibers = [FiberBlock()]
    for strong in scc:
        strong.check_input(graph)
        strong.classify_strong(graph)
        if strong.type==0:
            for node in strong.get_nodes():
                fibers[0].insert_node(node)
        elif strong.type==1:
            fibers.append(FiberBlock())
            for node in strong.get_nodes():
                fibers[-1].insert_node(node)
        elif strong.type == 2:
            node = strong.get_nodes()[0]
            fibers[0].insert_node(node)

    return fibers

def copy_class(copied_class):
    new_class = FiberBlock()
    for v in copied_class.get_nodes():
        new_class.insert_node(v)
    return new_class

def set_colors(graph, fibers):
    node_colors = graph.vp.node_colors
    
    for index, fclass in enumerate(fibers):
        for v in fclass.get_nodes():
            node_colors[v] = index

def set_ISCV(graph, ncolor):
    '''
        Define the ISCV for each node. At this step, each ISCV has size
        equal to 'ncolor' times the number of edge types (not considered yet).
    '''
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
            # 'Color_counter[k]' returns the number of inputs from color 'k'.
            iscv[v] += str(Colors_counter[k])

def split_fiber(fibernodes, fiber_iscv):
    #nfibers = len(Counter(fiber_iscv))
    iscv_order = np.argsort(fiber_iscv)

    new_list = []
    new_list.append(FiberBlock())
    current_iscv = fiber_iscv[iscv_order[0]]
    for index in iscv_order:
        if fiber_iscv[index] == current_iscv:
            new_list[-1].insert_node(fibernodes[index])
        else:
            current_iscv = fiber_iscv[index]
            new_list.append(FiberBlock())
            new_list[-1].insert_node(fibernodes[index])
    return new_list

def split_fiberf(class_index, fiber, fibernodes, fiber_iscv):
    '''
        Given the class that will be splitted, we delete from it
        the nodes that have different ISCV than the first element.

        This, if a class is splitted in N classes, then we will
        create N-1 new classes, and use the original one to store
        the class of the first element (in the ordered iscv).
    '''
    iscv_order = np.argsort(fiber_iscv)

    new_list = []
    new_list.append(FiberBlock())
    current_iscv = fiber_iscv[iscv_order[0]]
    for index in iscv_order:
        if fiber_iscv[index]!=current_iscv:
            current_iscv = fiber_iscv[index]
            new_list.append(FiberBlock())
        new_list[-1].insert_node(fibernodes[index])

    to_be_deleted = []
    for index, new_class in enumerate(new_list):
        if index==0: continue
        for v in new_class.get_nodes():
            to_be_deleted.append(v)
    fiber.delete_nodes(to_be_deleted)
    
    return new_list[1:]


def print_colors(graph):
    node_colors = graph.vp.node_colors
    for v in graph.get_vertices():
        print(v, node_colors[v])

def number_colors(graph):
    node_colors = list(graph.vp.node_colors.a)
    print(len(set(node_colors)))
