'''
    Based on the methods presented on the draft, we collect the
    runtime measurements for each fibration algorithm and we store
    on data files.
'''
import numpy as np
import graph_tool.all as gt
from fibration_mains import *

### From the "graph tool library" we create a random directed network ###
def deg_sampler():
    return (np.random.randint(1,16), np.random.randint(1,16))

N = 5000
g = gt.random_graph(N, deg_sampler, model="erdos", self_loops=True, random=True)
#########################################################################

# -> Minimal balanced coloring algorithm <- #
MBColoring(g)

# -> Fast fibration partitioning algorithm <- #
#regulation = g.new_edge_property('int')
#for n in g.edges():
#    regulation[n] = 0
#g.edge_properties['regulation'] = regulation

#FFPartitioning(g)


#pos_random = gt.random_layout(g)
#gt.graph_draw(g, pos=pos_random , output="../Figures/fibers/random_test.pdf")

