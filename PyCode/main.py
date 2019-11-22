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

g = buildGraph(edgefile)

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
while bqueue:
	print(len(bqueue))
	refinement_set = bqueue.popleft()
	INPUT_SPLIT(partition, refinement_set, g, bqueue)

print(len(partition))

count = 0
for eachblock in partition:
	eachblock.show_nodes()
	size = eachblock.get_number_nodes()
	if size>1: count+=1

print(count)


