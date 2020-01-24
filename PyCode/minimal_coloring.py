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
from minimalcoloringf import *
from collections import Counter

identifier = sys.argv[1]
edgefile = "../Data/"+identifier+"edgelist.dat"

g = buildGraph(edgefile)
N = g.get_vertices().shape[0]

#####################################################################
############### MINIMAL BALANCED COLORING ALGORITHM #################
#####################################################################

# Properties of the nodes: colors and ISCV.
node_colors = g.new_vertex_property('int')
iscv = g.new_vertex_property('string')
g.vertex_properties['node_colors'] = node_colors
g.vertex_properties['iscv'] = iscv

#### INITIALIZATION: Every node has the same color. ####
ncolor_before = -1
ncolor_after = set_SOLITAIRE(g) # Initial number of colors.
n_unique = ncolor_after         # Number of unique iscv.

#solitaires = get_SOLITAIRE(g)   # Solitaire nodes.  
########################################################

set_ISCV(g, ncolor_after)
print(list(g.vp.iscv))

#### Refinement loop ####
while ncolor_after!=ncolor_before:
    iscv_list = list(g.vp.iscv)
    iscv_count = Counter(iscv_list)
    n_unique = len(iscv_count)

    ncolor_before = ncolor_after
    ncolor_after = n_unique
    
    # Gets the unique ISCV and assign to them different colors.
    unique_iscv = list(set(iscv_count.elements()))
    colors_iscv = np.arange(0, n_unique, 1)

    # Each node receives a color label according to its ISCV label.
    for node in g.get_vertices():
        g.vp.node_colors[node] = colors_iscv[unique_iscv.index(g.vp.iscv[node])]
    for sol in solitaires[1:]:
        g.vp.node_colors[node] = n_unique
        n_unique+=1

    set_ISCV(g, ncolor_after) # putting zero strings "0000..." in the same class.

#number_colors(g)




    














