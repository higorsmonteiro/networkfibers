import sys
import numpy as np
from utils import *
from fiber import *
from FFPf import *
#from fibrationf import *
import graph_tool.all as gt
from collections import deque, defaultdict   # lists as queues are very slow.

identifier = sys.argv[1]
flagname = sys.argv[2]

edgefile = "../Data/"+identifier+"edgelist.dat"
nodename = "../Data/"+identifier+"nameID.dat"

if flagname=='-y':
    g, n_regulation = defineGraph(edgefile, nodenamefile=nodename)
else:
    g, n_regulation = defineGraph(edgefile)

N = g.get_vertices().shape[0]

############# COARSEST REFINEMENT PARTITIONING ALGORITHM ##############
''' 
    Define the initial partition as one block containing all operating 
    nodes. Operating nodes are nodes that receive information at least
    from itself, which not include nodes that do not receive any 
    information: the solitaire.
'''

bqueue = deque([])
partition = Initialization(g, bqueue)

# Until the queue is empty, we procedure the splitting process.
while bqueue:
    pivot_set = bqueue.popleft()
    #refinement_set.show_nodes()
    input_splitf(partition, pivot_set, g, n_regulation, bqueue)
    #INPUT_SPLIT(partition, pivot_set, g, bqueue)

PrintFibers(partition, g)
### Check input-set stability with respect to all fibers.
#regulation = g.edge_properties['regulation'].a
#first_fiber = partition[0]
#for index, block in enumerate(partition):
#    block.index = index
#    b = first_fiber.input_stability(g, block, regulation)
#    if b==(-1):
#        print("Not stable")
#        break
#########################################################
#PrintFibers(partition, g)

#fibercolors = g.new_vertex_property('float')
#fiberindex = g.new_vertex_property('int')
#
#for block in partition:
#    nodes = block.get_nodes()
#    for node in nodes: fiberindex[node] = block.index
#
#for block in solitaire:
#    nodes = block.get_nodes()
#    for node in nodes: fiberindex[node] = -1
#
#f = fiberindex.a
#print(f)
#
#vfilt = g.new_vertex_property('bool')
#vfilt.a[f==8] = True
#print(vfilt)
#
#f = gt.GraphView(g, vfilt=vfilt)
#
#pos_random = gt.random_layout(f)
#pos = gt.arf_layout(f, max_iter=0, a=0.2)
#gt.graph_draw(f, pos=pos_random, vertex_text=g.vertex_index, vertex_font_size=26, output="../Figures/fibers/test.pdf")   


################ Draw all the fibers ################ 
# Defines a graph structure to be used for each fiber.


#####################################################



#cc = GetNumberFibers(partition)
#print(cc)
#print("strong")
#st = GetStrongCCompOfNode(g, 3)
#print(st)
#print("dist")
#PrintFibers(partition, g)
#n = GetNumberFibers(partition)
#print(n)

#print("Input:")
#counting = GetNumberFibers(partition)
#if flagname=="-y": PrintFibers(partition, g, name=True)
#else: PrintFibers(partition, g)

#print("Output:")
#counting_out = GetNumberFibers(partition_out)
#if flagname=="-y": PrintFibers(partition_out, g, name=True)
#else: PrintFibers(partition_out, g)

#print(counting)
#print(counting_out)
