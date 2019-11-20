from fiber import *
from utils import *
import numpy as np
from collections import deque

def UPGRADE_PARTITION(new_blocks, old_blocks, partition):
    for old in old_blocks:
        partition.remove(old)

    for new in new_blocks:
        partition.append(new)


def PUSH_ON_BLOCK(node, indexlist, splitted_part):
    for fblock in splitted_part:
        if indexlist.sort() == fblock.regtype.sort():
            fblock.insert_node(node)
            return
    
    newblock = FiberBlock()
    newblock.regtype = indexlist
    newblock.insert_node(node)
    splitted_part.append(newblock)

def SPLIT_BLOCK(fiberblock, pos_fromSet, neg_fromSet, dual_fromSet, subpart2):
    splitted = []
    fibernodes = fiberblock.get_nodes()
    for node in fibernodes:
        PUSH_ON_BLOCK(node, [pos_fromSet[node], neg_fromSet[node], dual_fromSet[node]], splitted)

    index = 0
    block_sizes = []
    for fblock in splitted:
        fblock.index = index
        block_sizes.append(fblock.get_number_nodes())
        index += 1

    maxindex = block_sizes.index(max(block_sizes))
    splitted[maxindex].index = -96

    for block in splitted: subpart2.append(block)


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
        SPLIT_BLOCK(fblock, pos_fromSet, neg_fromSet, dual_fromSet, subpart2)


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

    BLOCKS_PARTITIONING(subpart1, subpart2, pos_fromSet, neg_fromSet, dual_fromSet)

    if len(subpart2)>len(subpart1):
        UPGRADE_PARTITION(subpart2, subpart1, partition)

        for splitted in subpart2:
            if splitted.index!=(-96): bqueue.append(splitted)

    
    



