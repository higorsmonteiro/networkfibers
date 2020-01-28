import numpy as np
import graph_tool.all as gt
from fibration_mains import *

N = 100

def deg_sampler():
    return (np.random.randint(1,16), np.random.randint(1,16))

g = gt.random_graph(N, deg_sampler, model="erdos", self_loops=True, random=True)

#MBColoring(g)
regulation = g.new_edge_property('int')
for n in g.edges():
    regulation[n] = 0
g.edge_properties['regulation'] = regulation

FFPartitioning(g)


#pos_random = gt.random_layout(g)
#gt.graph_draw(g, pos=pos_random , output="../Figures/fibers/random_test.pdf")

