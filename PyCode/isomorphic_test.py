import sys
import numpy as np
from utils import *
from fiber import *
from fibrationf import *
import graph_tool.all as gt
from collections import deque   # lists as queues are very slow.

identifier1 = sys.argv[1]
identifier2 = sys.argv[2]

edgefile1 = "../Data/"+identifier1+"edgelist.dat"
nodename1 = "../Data/"+identifier1+"nameID.dat"

edgefile2 = "../Data/"+identifier2+"edgelist.dat"
nodename2 = "../Data/"+identifier2+"nameID.dat"

g1 = buildGraph(edgefile1)
g2 = buildGraph(edgefile2)

if g1.get_vertices().shape[0]!=g2.get_vertices().shape[0]: 
    print("Non-Iso") 
    quit()
##################################################################

############ FOR FIRST NETWORK ############
partition1 = []
solitaire1 = []
bqueue1 = deque([])
PREPROCESSING(g1, partition1, solitaire1, bqueue1)
ENQUEUE_BLOCKS(partition1, bqueue1)
ENQUEUE_BLOCKS(solitaire1, bqueue1)

partition_out1 = []
solitaire_out1 = []
bqueue_out1 = deque([])
OUTPUT_PREPROCESSING(g1, partition_out1, solitaire_out1, bqueue_out1)
ENQUEUE_BLOCKS(partition_out1, bqueue_out1)
ENQUEUE_BLOCKS(solitaire_out1, bqueue_out1)

while bqueue1:
	refinement_set = bqueue1.popleft()
	INPUT_SPLIT(partition1, refinement_set, g1, bqueue1)

while bqueue_out1:
	refinement_set = bqueue_out1.popleft()
	OUTPUT_SPLIT(partition_out1, refinement_set, g1, bqueue_out1)

print("NETWORK 1:")
print("Input:")
counting = GetNumberFibers(partition1)
#PrintFibers(partition1, g1)

print("Output:")
counting_out = GetNumberFibers(partition_out1)
#PrintFibers(partition_out1, g1)
print(counting, counting_out)
#######################################################################

############ FOR SECOND NETWORK ############
partition2 = []
solitaire2 = []
bqueue2 = deque([])
PREPROCESSING(g2, partition2, solitaire2, bqueue2)
ENQUEUE_BLOCKS(partition2, bqueue2)
ENQUEUE_BLOCKS(solitaire2, bqueue2)

partition_out2 = []
solitaire_out2 = []
bqueue_out2 = deque([])
OUTPUT_PREPROCESSING(g2, partition_out2, solitaire_out2, bqueue_out2)
ENQUEUE_BLOCKS(partition_out2, bqueue_out2)
ENQUEUE_BLOCKS(solitaire_out2, bqueue_out2)

while bqueue2:
	refinement_set = bqueue2.popleft()
	INPUT_SPLIT(partition2, refinement_set, g2, bqueue2)

while bqueue_out2:
	refinement_set = bqueue_out2.popleft()
	OUTPUT_SPLIT(partition_out2, refinement_set, g2, bqueue_out2)

print("NETWORK 2:")
print("Input:")
counting = GetNumberFibers(partition2)
#PrintFibers(partition2, g2)

print("Output:")
counting_out = GetNumberFibers(partition_out2)
#PrintFibers(partition_out2, g2)
print(counting, counting_out)

