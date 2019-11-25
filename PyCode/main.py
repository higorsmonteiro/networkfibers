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

if flagname=='-y':
    g = buildGraph(edgefile, nodenamefile=nodename)
else:
    g = buildGraph(edgefile)

############# COARSEST REFINEMENT PARTITIONING ALGORITHM ##############

''' 
    Define the initial partition as one block containing all operating 
    nodes. Operating nodes are nodes that receive information at least
    from itself, which not include nodes that do not receive any 
    information: the solitaire.
'''

partition = []
solitaire = []
bqueue = deque([])

PREPROCESSING(g, partition, solitaire, bqueue)
ENQUEUE_BLOCKS(partition, bqueue)
ENQUEUE_BLOCKS(solitaire, bqueue)

# Until the queue is empty, we procedure the splitting process.
while bqueue:
	refinement_set = bqueue.popleft()
	INPUT_SPLIT(partition, refinement_set, g, bqueue)

### Check input-set stability with respect to all fibers.
#regulation = g.edge_properties['regulation'].a
#first_fiber = partition[0]
#for block in partition:
#	b = first_fiber.input_stability(g, block, regulation)
#	print(b)
#########################################################

counting = GetNumberFibers(partition)
if flagname=="-y": PrintFibers(partition, g, name=True)
else: PrintFibers(partition, g)

print(counting)
