from fiber import *
from utils import *
import numpy as np
import graph_tool.all as gt
from collections import deque, defaultdict

#######################################################################
def Initialization(graph, bqueue):
    ''' 
        Separate each node in its corresponding SCC using
        the 'label_components' function from graph_tool
        library. After define the SCC we classify each one
        to correctly initialize the fibration algorithm.
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
        node_list = fibers_listing[scc_j]
        for n in node_list: scc[scc_j].insert_node(n)
    '''
        'scc[j]' is an object containing the information
        about the nodes inside the j-th SCC.
    '''
    partition = [FiberBlock()]
    autopivot = []
    ''' Defines if each SCC receives or not input
        from other components not itself. '''
    for strong in scc:
        strong.check_input(graph)
        strong.classify_strong(graph)
        if strong.type == 0:    # receive external input.
            for node in strong.get_nodes():
                partition[0].insert_node(node)
        elif strong.type == 1:  # SCC does not receive any external input.
            partition.append(FiberBlock())
            for node in strong.get_nodes():
                partition[-1].insert_node(node)
        elif strong.type == 2:  # does not receive external input, but it is an isolated autorregulated node.
            node = strong.get_nodes()[0]
            #for node in strong.get_nodes():
            autopivot.append(FiberBlock())
            autopivot[-1].insert_node(node)
            partition[0].insert_node(node)
        

    fiber_index = graph.vp.fiber_index
    for index, init_class in enumerate(partition): 
        bqueue.append(copy_class(init_class))
        for v in init_class.get_nodes(): fiber_index[v] = index
    for isolated in autopivot: bqueue.append(copy_class(isolated))

    return partition
#########################################################################

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

def enqueue_splitted(newclass_list, bqueue, graph):
    class_sizes = np.array([fclass.get_number_nodes() for fclass in newclass_list])
    argmax = np.argmax(class_sizes)

    for index, fclass in enumerate(newclass_list):
        if index==argmax: continue
        new_pivot = copy_class(fclass)
        bqueue.append(new_pivot)

def class_split(chi, R_class, g, graph, partition):
    '''
        'R_class' is an unstable matrix in which its columns
        are not all equal. To split correctly the equivalent class
        we must sort this matrix according each column sequence of
        it. This way, columns that are equal will gonna be side by
        side, which facilitates the creation of the new classes.
    '''
    N_class = chi.get_number_nodes()
    Z = chi.get_nodes()

    # Really fast splitting procedure.
    nclass_list = []
    check_uniquecol = defaultdict(lambda:-1)
    for n in range(N_class):
        cur_str = ''.join(str(num) for num in R_class[:,n])
        col_index = check_uniquecol[cur_str]
        if col_index==-1:
            nclass_list.append(FiberBlock())
            check_uniquecol[cur_str] = len(nclass_list)-1
            nclass_list[-1].insert_node(Z[n])
        else:
            nclass_list[col_index].insert_node(Z[n])
    ########################################################
    
    #str_input = [np.array2string(R_class[:,n], separator="")[1:-1] for n in range(N_class)]
    #ORDER = np.argsort(str_input)
    #nclass_list = []
    #nclass_list.append(FiberBlock())
    #current_str = str_input[ORDER[0]]
    #for index in ORDER:
    #    if str_input[index]!=current_str:
    #        current_str = str_input[index]
    #        nclass_list.append(FiberBlock())    
    #    nclass_list[-1].insert_node(Z[index])

    class_index = graph.vertex_properties['fiber_index'][Z[0]]

    to_be_deleted = []
    for index, new_class in enumerate(nclass_list):
        if index==0: continue # the nodes of the first class will remain in the original class. 
        partition.append(copy_class(new_class))
        for v in new_class.get_nodes():
            to_be_deleted.append(v)
            graph.vertex_properties['fiber_index'][v] = len(partition) - 1
    partition[class_index].delete_nodes(to_be_deleted)
    
    return nclass_list

def fast_partitioning(receiver_classes, eta, f, R, partition, n_edgetype, bqueue, graph):
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
            newclass_list = class_split(chi, R_class, g, graph, partition)
            enqueue_splitted(newclass_list, bqueue, graph)


def input_splitf(partition, pivot, graph, n_edgetype, bqueue):
    eta = pivot.sucessor_nodes(graph)
    regulation = graph.edge_properties['regulation'].a
    # Given the node number, 'f' gives its index in 'eta'.
    f = defaultdict(lambda:-1)
    for eta_index, sucessor in enumerate(eta):  f[sucessor] = eta_index

    # 'R' represents a matrix (n_edgetype, len(eta)).
    R = np.vstack([np.zeros(len(eta), int) for row in range(n_edgetype)])

    calc_R(R, graph, pivot, f, regulation)
    receiver_classes = get_possible_unstable_classes(graph, pivot, partition)
    fast_partitioning(receiver_classes, eta, f, R, partition, n_edgetype, bqueue, graph)
    
