import sys
import numpy as np
from utils import *
from fiber import *
from fibrationf import *
import graph_tool.all as gt
from collections import deque   # lists as queues are very slow.

identifier = sys.argv[1]
flagname = sys.argv[2]

edgefile = "../Data/"+identifier+"edgelist.dat"
nodename = "../Data/"+identifier+"nameID.dat"

g, regulation = buildGraph(edgefile)


############# COARSEST REFINEMENT PARTITIONING ALGORITHM ##############

''' 
    Define the initial partition as one block containing all operating 
    nodes. Operating nodes are nodes that receive information at least
    from itself, which not include nodes that do not receive any 
    information: the solitaire.
'''

partition = []
solitaire_part = []
bqueue = deque([])

PREPROCESSING(g, partition, solitaire_part)
ENQUEUE_BLOCKS(partition, bqueue)
ENQUEUE_BLOCKS(solitaire_part, bqueue)

# Until the queue is empty, we procedure the splitting process.
#while bqueue:
#    refinement_set = bqueue.popleft()
#    INPUT_SPLIT(partition, refinement_set, g, bqueue)

