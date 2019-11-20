from fiber import *
from utils import *
import numpy as np
from collections import deque

def PREPROCESSING(graph, partition, solitaire_part):
    init_block = FiberBlock()

    for node in graph.get_vertices():
        solitaire_bool = IDENTIFY_SOLITAIRE(graph, node)

        if solitaire_bool==0:
            block = FiberBlock()
            block.insert_node(node)
            solitaire_part.append(block)
        elif solitaire_bool==1:
            block = FiberBlock()
            block.insert_node(node)
            partition.append(block)
        else:
            init_block.insert_node(node)
    
    partition.append(init_block)

def ENQUEUE_BLOCKS(partition, block_queue):
    for fiberblock in partition:
        block_queue.append(fiberblock)

def GET_NONSTABLE_BLOCKS(partition, subpart, pos_fromSet, neg_fromSet, dual_fromSet):
    for fblock in partition:
        fibernodes = fblock.get_nodes()

        poslist = pos_fromSet[fibernodes]
        neglist = neg_fromSet[fibernodes]
        dualist = dual_fromSet[fibernodes]
        if np.count_nonzero(poslist==poslist[0]) != poslist.shape[0]:
            subpart.append(fblock)
        elif np.count_nonzero(neglist==neglist[0]) != neglist.shape[0]:
            subpart.append(fblock)
        elif np.count_nonzero(dualist==dualist[0]) != dualist.shape[0]:
            subpart.append(fblock)

def BLOCKS_PARTITIONING(subpart1, subpart2, pos_fromSet, neg_fromSet, dual_fromSet):
    for fblock in subpart1:
        fibernodes = fblock.get_nodes()
        poslist = pos_fromSet[fibernodes]
        neglist = neg_fromSet[fibernodes]
        dualist = dual_fromSet[fibernodes]



def INPUT_SPLIT(partition, refinement_set, graph, bqueue):
    subpart1 = []
    subpart2 = []

    pos_fromSet = np.array(graph.get_vertices().shape[0], int)
    neg_fromSet = np.array(graph.get_vertices().shape[0], int)
    dual_fromSet = np.array(graph.get_vertices().shape[0], int)

    for node in graph.get_vertices():
        edgefromSet(pos_fromSet, graph, node, refinement_set, 0)
        edgefromSet(neg_fromSet, graph, node, refinement_set, 1)
        edgefromSet(dual_fromSet, graph, node, refinement_set, 2)

    GET_NONSTABLE_BLOCKS(partition, subpart1, pos_fromSet, neg_fromSet, dual_fromSet)

    
    



