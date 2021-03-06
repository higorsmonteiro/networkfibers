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

def set_ISCV(graph, ncolor, num_edgetype):
    '''
        Define the ISCV for each node in the graph. At this step, 
        each ISCV has size equal to 'ncolor' times 'num_edgetype'.

        The 'iscv' internal property of each node is a string type,
        where this string is the sequence given by an array 'iscv_v'.
        The construction of 'iscv_v' is defined by an array of size
        'ncolor' x 'num_edgetype', and considering that a node 'v'
        receives an edge of color 'n' and with type index 'ne', the 
        entry of index 'ncolor*ne + n' is increased by an unit.

        Finally 'iscv' receives the string format of 'iscv_v'.
    '''
    # internal graph properties
    iscv = graph.vp.iscv
    fiber_index = graph.vp.fiber_index
    regulation = graph.ep.regulation.a

    # define the iscv for each node.
    for v in graph.get_vertices():
        type_edge = []              # type index of each incoming edge.
        input_colors = []           # color of each incoming edge.
        # each row in 'in_edges' corresponds to ['source', 'target', 'edge_index']
        in_edges = graph.get_in_edges(v, [graph.edge_index])
        for in_neigh in in_edges:
            input_colors.append(fiber_index[in_neigh[0]])
            type_edge.append(regulation[in_neigh[2]])

        iscv_v = np.zeros(ncolor*num_edgetype, int) # ISCV(v)
        for ind_color, color in enumerate(input_colors):
            etype = type_edge[ind_color]
            iscv_v[ncolor*etype + color] += 1
        iscv[v] = ''.join(str(num) for num in iscv_v)
        #iscv[v] = np.array2string(iscv_v, separator="")[1:-1]
        

def split_fiberf(class_index, fiber, fibernodes, fiber_iscv):
    '''
        'fiber_iscv' is a list of strings, where the strings 
        represents the ISCVs of each node of the given class.
        From this list, this function determines the correct
        splitting.
        
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
#########################################################################

###### FUNCTIONS USED FOR CHECKING THE CORRECTNESS OF THE RESULTS ####### 
def get_possible_unstable_classes(graph, pivot, partition):
    ''' 
        We define the list A of classes that receives information from
        'pivot'. To guarantee that the fibration algorithm runs in loglinear
        time, this procedure should take advantage of the 'fiber_index'
        vertex property in the network to avoid a sweep over all the nodes
        in it. 'defineGraph' and 'fast_gnp_random' functions in 'utils.py' 
        guarantees this vertex property.

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
    ''' 
        For each class 'chi' in 'receiver_classes' that receives
        information from pivot, we define a 'R_class' matrix with
        size (n_edgetype, N_class) where N_class is the number of 
        nodes in 'chi'. We fill 'R_class' according to matrix 'R'.
   '''
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
            return chi
        return None

def check_sucessor_stability(pivot, fibers, graph, n_edgetype):
    '''
        Given two classes 'fclass' and 'pivot', this function returns
        True if fclass is input-set stable with respect to pivot.
    '''
    eta  = pivot.sucessor_nodes(graph)
    regulation = graph.ep.regulation.a

    f = defaultdict(lambda:-1)
    for eta_index, sucessor in enumerate(eta): f[sucessor] = eta_index

    R = np.vstack([np.zeros(len(eta), int) for row in range(n_edgetype)])
    calc_R(R, graph, pivot, f, regulation)
    receiver_classes = get_possible_unstable_classes(graph, pivot, fibers)
    chi = fast_checking(receiver_classes, eta, f, R, fibers, n_edgetype, graph)
    if chi!=None:
        print(-1, pivot.get_nodes(), chi.get_nodes())

################################################################################

def print_colors(graph):
    fiber_index = graph.vp.fiber_index
    for v in graph.get_vertices():
        print(v, fiber_index[v])

def number_colors(graph):
    fiber_index = list(graph.vp.fiber_index.a)
    print(len(set(fiber_index)))

def copy_class(copied_class):
    new_class = FiberBlock()
    for v in copied_class.get_nodes():
        new_class.insert_node(v)
    return new_class
