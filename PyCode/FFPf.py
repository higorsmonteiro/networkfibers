from fiber import *
from utils import *
import numpy as np
import graph_tool.all as gt
from collections import deque, defaultdict

#######################################################################
def Initialization(graph, nodes_obj, bqueue):
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
        elif strong.type == 2:  # does not receive external input, but it is an isolated autorregulated node.
            for node in strong.get_nodes():
                fiber = FiberBlock()
                fiber.insert_node(node)
                partition[0].insert_node(node)
                autopivot.append(fiber)
                #bqueue.append(fiber)
        elif strong.type == 1:  # SCC does not receive any external input.
            partition.append(FiberBlock())
            for node in strong.get_nodes():
                partition[-1].insert_node(node)

    for index, fiber in enumerate(partition): 
        bqueue.append(fiber)
        for v in fiber.get_nodes(): nodes_obj[v].set_class_index(index)
    
    for isolated in autopivot: bqueue.append(isolated)
    return partition
#########################################################################


def upgrade_partition(new_classes, old_classes, partition):
    for old in old_classes: partition.remove(old)
    for new in new_classes: partition.append(new)

#def efficient_get_unstable_classes(partition, unstable_list, regulation_list, pivotnode_to_index):


def get_unstable_classes(partition, unstable_list, regulation_list, pivotnode_to_index):
    ''' it is not linerar '''
    for fclass in partition:
        class_size = fclass.get_number_nodes()
        nodelist = fclass.get_nodes()
        if class_size==0: continue

        first_node = nodelist[0]
        for type_fromSet in regulation_list:
            for k in range(1, class_size):
                first_node_index = pivotnode_to_index[first_node]    
                correct_index = pivotnode_to_index[nodelist[k]]
                if first_node_index == -1 and correct_index != -1:
                    if type_fromSet[correct_index]!=0:
                        unstable_list.append(fclass)
                        break
                elif correct_index == -1 and first_node_index != -1:
                    if type_fromSet[first_node_index]!=0:
                        unstable_list.append(fclass)
                        break
                elif correct_index == -1 and first_node_index == -1:
                    pass
                elif type_fromSet[correct_index]!=type_fromSet[first_node_index]:
                    unstable_list.append(fclass)
                    break

def push_to_class(node, indexlist, splitted_part):
    ''' 'indexlist' represents the sequence of number of
        connections received from the pivot set for each 
        edge type   ''' 
    #print(indexlist)
    seqlist = np.array(indexlist)
    for fclass in splitted_part:
        class_types = np.array(fclass.regtype)
        if np.array_equal(seqlist, class_types):
            fclass.insert_node(node)
            return

    newclass = FiberBlock()
    newclass.regtype = indexlist
    newclass.insert_node(node)
    splitted_part.append(newclass)

def classes_partitioning(unstable_classes, splitted_classes, regulation_list, node_to_index):
    
    for fclass in unstable_classes:
        splitted = []
        nodelist = fclass.get_nodes()
        for node in nodelist:
            nodeindex = node_to_index[node]
            if nodeindex==-1: regulation_node = [0 for type_arr in regulation_list]
            else: regulation_node = [type_arr[nodeindex] for type_arr in regulation_list]

            #print(node, nodeindex, len(regulation_list[:]))
            push_to_class(node, regulation_node, splitted)

        #print(len(splitted))
        index = 0
        classes_size = []
        for sclass in splitted:
            sclass.index = index
            classes_size.append(sclass.get_number_nodes())
            splitted_classes.append(sclass)
            index += 1
        
        maxindex = classes_size.index(max(classes_size))
        splitted[maxindex].index = -96


def possible_unstable_c(partition, pivot_sucessors):
    maybe_unstables = []
    for fclass in partition:
        class_nodes = fclass.get_nodes()
        bool_list = [j for j in pivot_sucessors if j in class_nodes]
        if len(bool_list)>0:
            maybe_unstables.append(fclass)
    return maybe_unstables

def input_splitf(partition, pivot, graph, n_edgetype, bqueue):
    ''' The splitting process is divided in the following steps:

        1.  We get all the current classes that are unstable with
            respect to 'pivot'. We do that by getting all the outgoing
            neighbors of the pivot nodes and their classes. From these
            classes, we select only the ones that are unstable.

        2.  

    '''
    unstable_classes = []
    splitted_classes = []
    N = graph.get_vertices().shape[0]
    
    pivot_sucessors = pivot.sucessor_nodes(graph)
    # Given the node number, 'node_to_index' gives its index in 'pivot_sucessors'.
    node_to_index = defaultdict(lambda:-1)
    for index, sucessor in enumerate(pivot_sucessors):
        node_to_index[sucessor] = index

    # 'regulation_list' represents a matrix (n_edgetype, len(pivot_nodes)).
    regulation_list = []   
    for j in range(n_edgetype):
        regulation_list.append(np.zeros(len(pivot_sucessors), int))

    ''' All nodes that receives information from 'pivot'
        receives, for each edge type, the number of incoming
        links received from 'pivot'. '''
        # Type of each edge: 0,1,2,...,# edge types - 1.
    regulation = graph.edge_properties['regulation'].a
    edgefromSet_optimal(regulation_list, graph, pivot, node_to_index, regulation)
    possibles = possible_unstable_c(partition, pivot_sucessors)

    get_unstable_classes(possibles, unstable_classes, regulation_list, node_to_index)
    classes_partitioning(unstable_classes, splitted_classes, regulation_list, node_to_index)
    
    if len(splitted_classes)>len(unstable_classes):
        upgrade_partition(splitted_classes, unstable_classes, partition)
        for splitted in splitted_classes:
            if splitted.index!=(-96): bqueue.append(splitted)



