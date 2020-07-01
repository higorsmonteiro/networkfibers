import sys
import os
import numpy as np
import pandas as pd
import graph_tool.all as gt
from collections import defaultdict, deque

from lib.fiber import FiberBlock, StrongComponent

class FFIBRATION:
    def __init__(self):
        self.g = None
        self.initial_partition = None
        self.queue = None

    def get_network(self):
        return self.g
        
    def generate_network(self, data_path, fname, sep=' ', def_dict=True, weight=False):
        
        # Read the data with an appropriate format.
        table = pd.read_csv(os.path.join(data_path, fname), sep=sep, header=None)
        m = table.shape[0]
        ncols = table.shape[1]
        
        src_nodes = table[0] # Source node labels.
        tgt_nodes = table[1] # Target node labels.
        # If there is an extra column for categorical weights.
        if ncols>2 and weight:
            link_info = table[2]
            weight_unique = link_info.unique()
            w_label = defaultdict(lambda:-1)
            
            n = 0
            for weight_name in weight_unique:
                w_label[weight_name] = n
                n += 1
        else:
            # One single category.
            table['weight_label'] = table[2].apply(lambda x:0)
            
        # Initialize the network
        G = gt.Graph(directed=True)
        
        wlink_values = G.new_ep('int')

        node_unique = np.unique(np.concatenate((src_nodes, tgt_nodes)))
        N = node_unique.shape[0]

        # create dictionary with the labels for each node.
        n = 0
        node_index = defaultdict(lambda:-1)
        for v_name in node_unique:
            if node_index[v_name]==-1:
                node_index[v_name] = n
                n += 1

        table['src_label'] = table[0].apply(lambda x: node_index[x])
        table['tgt_label'] = table[1].apply(lambda x: node_index[x])
        sources = np.array(table['src_label'])
        targets = np.array(table['tgt_label'])
        if ncols>2 and weight:
            table['weight_label'] = table[2].apply(lambda x: w_label[x])

        
        weights = np.array(table['weight_label'])
        fmt_edgelist = np.vstack((sources, targets, weights)).T
        G.add_edge_list(fmt_edgelist, eprops=[wlink_values])

        wlink_labels = G.new_ep('string')
        for k, e in enumerate(G.edges()): wlink_labels[e] = table[2].iloc[k]
        G.edge_properties['link_name'] = wlink_labels
        G.edge_properties['link_type_index'] = wlink_values

        fiber_index = G.new_vp('int')
        fiber_index.a = np.zeros(N, int)
        node_labels = G.new_vp('string')
        for key in node_index.keys():
            node_labels[node_index[key]] = key
        G.vertex_properties['node_name'] = node_labels
        G.vertex_properties['fiber_index'] = fiber_index
        self.g = G.copy()

    def initial_setup(self):
        bqueue = deque([])
        partition = initialization(self.g, bqueue)

        self.initial_partition = partition
        self.queue = bqueue


## ------------------- FUNCTIONS --------------------- ##

def initialization(graph, bqueue):
    ''' 
        Separate each node in its corresponding SCC using
        the 'label_components' function from graph_tool
        library. After defining the SCC we classify each one
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
            autopivot.append(FiberBlock())
            autopivot[-1].insert_node(node)
            partition[0].insert_node(node)

    fiber_index = graph.vp.fiber_index
    for index, init_class in enumerate(partition): 
        bqueue.append(copy_class(init_class))
        for v in init_class.get_nodes(): fiber_index[v] = index
    for isolated in autopivot: bqueue.append(copy_class(isolated))

    return partition

def copy_class(copied_class):
    new_class = FiberBlock()
    for v in copied_class.get_nodes():
        new_class.insert_node(v)
    return new_class


        




    