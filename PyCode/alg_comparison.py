'''
    Based on the methods presented on the draft, we collect the
    runtime measurements for each fibration algorithm and we store
    on data files.
'''
import sys
import numpy as np
import graph_tool.all as gt
from fibration_mains import *

### From the "graph tool library" we create a random directed network ###
def deg_sampler():
    return (np.random.randint(1,10), np.random.randint(1,10))

N = int(sys.argv[1])
mode = sys.argv[2]

g = gt.random_graph(N, deg_sampler, model="erdos", self_loops=True, random=True)
regulation = g.new_edge_property('int')
for n in g.edges():
    regulation[n] = 0
g.edge_properties['regulation'] = regulation
print(g.get_edges().shape[0])
#########################################################################

if mode=="-mbc":
    # -> Minimal balanced coloring algorithm <- #
    MBColoring(g)
elif mode=="-fpp":
    # -> Fast fibration partitioning algorithm <- #
    FFPartitioning(g)
else:
    print("Wrong mode flag.")


