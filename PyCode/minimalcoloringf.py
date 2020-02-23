import numpy as np
from utils import *
from fiber import *
import graph_tool.all as gt
from collections import Counter, defaultdict

def Initialization(graph):
    ''' 
        Separate each node in its corresponding SCC and WCC 
        using the 'label_components' function from graph_tool
        library.
    '''
    #weak_components = gt.label_components(g, directed=False)
    #graph.vertex_properties['weak_components'] = weak_components
    
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

def get_possible_unstable_classes(graph, pivot, partition):
    ''' 
        We define the list A of classes that receives information from
        'pivot'. To guarantee that the fibration algorithm runs in loglinear
        time, this procedure should take advantage of the 'fiber_index'
        vertex property in the network to avoid a sweep over all the nodes
        in it. 'defineGraph' function in 'utils.py' guarantees this vertex
        property.

        By getting the fiber indexes of all outcoming neighbors of the pivot
        set nodes the 'receiver_classes' stores the indexes of the possible
        unstable classes. Then, 'classes' lists the corresponding objects. 

    '''
    fiber_index = graph.vp.fiber_index

    receiver_classes = []
    p_sucessors = pivot.sucessor_nodes(graph)
    for w in p_sucessors: receiver_classes.append(fiber_index[w])
    receiver_classes = set(receiver_classes)
    
    classes = [partition[f_index] for f_index in receiver_classes]
    return classes

def fast_checking(receiver_classes, eta, f, R, partition, n_edgetype, graph):
    ''' For each class 'chi' in 'receiver_class' that receives
        information from pivot, we define a 'R_class' matrix with
        size (n_edgetype, N_class) where N_class is the number of 
        nodes in 'chi'. We fill 'R_class' accordigly to matrix 'R'.   '''

    for chi in receiver_classes:
        chi_nodes = chi.get_nodes()
        N_class = chi.get_number_nodes()
        
        g = defaultdict(lambda:-1)
        for n in range(N_class): g[chi_nodes[n]] = n
        R_class = np.vstack([np.zeros(N_class, int) for j in range(n_edgetype)])
        
        # fill 'R_class' according 'R'. Efficient.
        for v in eta:
            if g[v]!=-1:
                for m in range(n_edgetype):
                    R_class[m,g[v]] = R[m,f[v]]
        
        if is_unstable(R_class):
            print(-1)
        #print(1)

def set_colors(graph, fiber_list):
    '''
        Receives a list of fibers, and for each node
        associate it with the index number of its fiber
        location on the list.
    '''
    fiber_index = graph.vp.fiber_index
    for index, fclass in enumerate(fiber_list):
        for v in fclass.get_nodes():
            fiber_index[v] = index

def set_ISCV(graph, ncolor):
    '''
        Define the ISCV for each node. At this step, each ISCV has size
        equal to 'ncolor' times the number of edge types (not considered yet).
    '''
    # intrinsic properties
    iscv = graph.vp.iscv
    fiber_index = graph.vp.fiber_index

    # define the iscv for each node.
    for v in graph.get_vertices():
        iscv[v] = ""
        in_neighbors = graph.get_in_neighbors(v)
        input_colors = []
        for neigh in in_neighbors:
            input_colors.append(fiber_index[neigh])
        # Counts how many inputs the 'node' receives from each color.
        Colors_counter = Counter(input_colors)
        # 'Color_counter[k]' returns the number of inputs from color 'k'.
        for k in range(ncolor): iscv[v] += str(Colors_counter[k])  

#def split_fiber(fibernodes, fiber_iscv):
#    iscv_order = np.argsort(fiber_iscv)
#
#    new_list = []
#    new_list.append(FiberBlock())
#    current_iscv = fiber_iscv[iscv_order[0]]
#    for index in iscv_order:
#        if fiber_iscv[index] == current_iscv:
#            new_list[-1].insert_node(fibernodes[index])
#        else:
#            current_iscv = fiber_iscv[index]
#            new_list.append(FiberBlock())
#            new_list[-1].insert_node(fibernodes[index])
#    return new_list

def split_fiberf(class_index, fiber, fibernodes, fiber_iscv):
    '''
        Given the class that will be splitted, we delete from it
        the nodes that have different ISCV than the first element.

        This, if a class is splitted in N classes, then we will
        create N-1 new classes, and use the original one to store
        the class of the first element (in the ordered iscv).
    '''
    iscv_order = np.argsort(fiber_iscv)

    new_list = [FiberBlock()]
    current_iscv = fiber_iscv[iscv_order[0]]
    for index in iscv_order:
        if fiber_iscv[index]!=current_iscv:
            current_iscv = fiber_iscv[index]
            new_list.append(FiberBlock())
        new_list[-1].insert_node(fibernodes[index])

    to_be_deleted = []
    for index, new_class in enumerate(new_list):
        if index==0: continue
        for v in new_class.get_nodes(): to_be_deleted.append(v)
    fiber.delete_nodes(to_be_deleted)
    return new_list[1:]


def print_colors(graph):
    fiber_index = graph.vp.fiber_index
    for v in graph.get_vertices():
        print(v, fiber_index[v])

def number_colors(graph):
    fiber_index = list(graph.vp.fiber_index.a)
    print(len(set(fiber_index)))
