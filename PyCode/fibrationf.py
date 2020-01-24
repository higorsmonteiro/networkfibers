from fiber import *
from utils import *
import numpy as np
from collections import deque

def UPGRADE_PARTITION(new_blocks, old_blocks, partition):
    for old in old_blocks:
        partition.remove(old)

    for new in new_blocks:
        partition.append(new)


def PREPROCESSING(graph, partition, solitaire_part, bqueue):
    '''
        'partition' will consist of all operating nodes,
        including the ones that only receive information
        from themselves. 'solitaire_part' will consist of
        all nodes that do not receive any information. The
        nodes that receive information only from themselves
        will be put in queue as single blocks.
    '''
    init_block = FiberBlock()
    all_nodes = graph.get_vertices()

    for node in all_nodes:
        solitaire_bool = IDENTIFY_SOLITAIRE(graph, node)
        
        if solitaire_bool==0: # Full solitaire.
            block = FiberBlock()
            block.insert_node(node)
            solitaire_part.append(block)
        elif solitaire_bool==1: # Receives information only from itself.
            block = FiberBlock()
            block.insert_node(node)
            bqueue.append(block)
            init_block.insert_node(node)
        else:   # Otherwise.
            init_block.insert_node(node)
    
    partition.append(init_block)

def ENQUEUE_BLOCKS(partition, block_queue):
    for fiberblock in partition:
        block_queue.append(fiberblock)

def GET_NONSTABLE_BLOCKS(partition, subpart, pos_fromSet, neg_fromSet, dual_fromSet):
    '''
        Get all the fiber blocks that are input-tree unstable with respect to a
        refinement set. To this, the function checks if, for all the edge types
        and for each fiber, the number of the edges received from the refinement set
        is the same.
    '''
    for fblock in partition:
        n = fblock.get_number_nodes()
        fibernodes = fblock.get_nodes()
        if n==0: continue

        first_node = fibernodes[0]
        for k in range(1, n):
            if pos_fromSet[first_node]!=pos_fromSet[fibernodes[k]]:
                subpart.append(fblock)
                break
            if neg_fromSet[first_node]!=neg_fromSet[fibernodes[k]]:
                subpart.append(fblock)
                break
            if dual_fromSet[first_node]!=dual_fromSet[fibernodes[k]]:
                subpart.append(fblock)
                break


def PUSH_ON_BLOCK(node, indexlist, splitted_part):
    node_types = np.array(indexlist)

    for fblock in splitted_part:
        block_types = np.array(fblock.regtype)
        if np.array_equal(node_types, block_types):
            fblock.insert_node(node)
            return
    
    newblock = FiberBlock()
    newblock.regtype = indexlist
    newblock.insert_node(node)
    splitted_part.append(newblock)

def BLOCKS_PARTITIONING(subpart1, subpart2, pos_fromSet, neg_fromSet, dual_fromSet):
    
    for fblock in subpart1:
        splitted = []
        fibernodes = fblock.get_nodes()
        
        # splitting process: 'splitted' will hold the splitted blocks from subpart1.
        for node in fibernodes:
            PUSH_ON_BLOCK(node, [pos_fromSet[node], neg_fromSet[node], dual_fromSet[node]], splitted)

        # indexation and size of each splitted block.
        index = 0
        blocks_size = []
        for block in splitted:
            block.index = index
            blocks_size.append(block.get_number_nodes())
            subpart2.append(block)
            index += 1

        maxindex = blocks_size.index(max(blocks_size))
        splitted[maxindex].index = -96


def INPUT_SPLIT(partition, refinement_set, graph, bqueue):
    subpart1 = []
    subpart2 = []

    regulation = graph.edge_properties['regulation'].a
    pos_fromSet = np.zeros(graph.get_vertices().shape[0], int)
    neg_fromSet = np.zeros(graph.get_vertices().shape[0], int)
    dual_fromSet = np.zeros(graph.get_vertices().shape[0], int)

    edgefromSet([pos_fromSet, neg_fromSet, dual_fromSet], graph, refinement_set, regulation)

    GET_NONSTABLE_BLOCKS(partition, subpart1, pos_fromSet, neg_fromSet, dual_fromSet)
    BLOCKS_PARTITIONING(subpart1, subpart2, pos_fromSet, neg_fromSet, dual_fromSet)

    if len(subpart2)>len(subpart1):
        UPGRADE_PARTITION(subpart2, subpart1, partition)

        for splitted in subpart2:
            if splitted.index!=(-96): bqueue.append(splitted)

############################################################################

def OUTPUT_PREPROCESSING(graph, partition, solitaire_part, bqueue):
    init_block = FiberBlock()
    all_nodes = graph.get_vertices()

    for node in all_nodes:
        solitaire_bool = OUT_IDENTIFY_SOLITAIRE(graph, node)
        
        if solitaire_bool==0: # Full solitaire.
            block = FiberBlock()
            block.insert_node(node)
            solitaire_part.append(block)
        elif solitaire_bool==1: # Sends information only to itself.
            block = FiberBlock()
            block.insert_node(node)
            bqueue.append(block)
            init_block.insert_node(node)
        else:   # Otherwise.
            init_block.insert_node(node)
    
    partition.append(init_block)

def OUTPUT_SPLIT(partition, refinement_set, graph, bqueue):
    subpart1 = []
    subpart2 = []

    regulation = graph.edge_properties['regulation'].a
    pos_toSet = np.zeros(graph.get_vertices().shape[0], int)
    neg_toSet = np.zeros(graph.get_vertices().shape[0], int)
    dual_toSet = np.zeros(graph.get_vertices().shape[0], int)

    edgetoSet([pos_toSet, neg_toSet, dual_toSet], graph, refinement_set, regulation)

    GET_NONSTABLE_BLOCKS(partition, subpart1, pos_toSet, neg_toSet, dual_toSet)
    BLOCKS_PARTITIONING(subpart1, subpart2, pos_toSet, neg_toSet, dual_toSet)

    if len(subpart2)>len(subpart1):
        UPGRADE_PARTITION(subpart2, subpart1, partition)

        for splitted in subpart2:
            if splitted.index!=(-96): bqueue.append(splitted)    



