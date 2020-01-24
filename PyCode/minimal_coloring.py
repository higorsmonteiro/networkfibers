'''
    Code for fiber identification using the minimal balanced coloring
    algorithm. We must pass as argument the identifier string for
    the file containing the edgelist of the network. The identifier
    does not contain the extension of the file, which is assumed to be
    .dat. This way, only the name is passed.  
'''

import sys
import numpy as np
from utils import *
import graph_tool.all as gt
from collections import Counter

identifier = sys.argv[1]
edgefile = "../Data/"+identifier+"edgelist.dat"

g = buildGraph(edgefile)
N = g.get_vertices().shape[0]

############### MINIMAL BALANCED COLORING ALGORITHM #################
ncolor = 1             # Number of colors.
n_unique = 1           # Number of unique iscv.

node_colors = g.new_vertex_property('int')
iscv = g.new_vertex_property('string')

#### INITIALIZATION: Every node has the same color. ####
for node in g.get_vertices():
    node_colors[node] = 0
    iscv[node] = ""

# Define the ISCV for each node. At this step, each ISCV has size 'ncolor'.
for node in g.get_vertices():
    # verifies all incoming neighbors
    in_neighbors = g.get_in_neighbors(node)
    input_colors = []
    for neigh in in_neighbors:
        input_colors.append(node_colors[neigh])
    # Counts how many inputs the 'node' receives from each color.
    Colors_counter = Counter(input_colors)
    for k in range(ncolor):
        # 'Color_counter[k] returns the number of inputs from color 'k'.
        iscv[node] += str(Colors_counter[k])
##########################################################
    














