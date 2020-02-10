'''
    Based on the methods presented on the draft, we collect the
    runtime measurements for each fibration algorithm and we store
    on data files.
'''
import sys
import timeit
import numpy as np
import graph_tool.all as gt

def fast_gnp_erdos(n, p, seed=None, gdirected=False):
    ''' fast procedure to generate an Erdos-Renyi network 
        according the G(n,p) model. '''
    if seed!=None:  np.seed(seed)

    G = gt.Graph(directed=gdirected)
    for nn in range(n): G.add_vertex()

    if p<=0 or p >= 1:  return None
    
    w = -1
    lp = np.log(1.0-p)
    if gdirected:   # directed network (self-loops allowed).
        v = 0
        while v < n:
            lr = np.log(1.0 - np.random.uniform())
            w += (1 + int(lr/lp))
            while v < n <= w:
                w -= n
                v += 1
            if v < n:
                G.add_edge(v,w)
    else:   # undirected network.
        v = 1
        while v < n:
            lr = np.log(1.0 - np.random.random())
            w += 1 + int(lr/lp)
            while w >= v and v < n:
                w -= v
                v += 1
            if v < n:
                G.add_edge(v, w)
    
    return G
#########################################################################

N = int(sys.argv[1])
p = float(sys.argv[2])
mode = sys.argv[3]
n_repeat = 2

######## FAST FIBRATION PARTITIONING ########
def FFP_time():
    SETUP_CODE = '''
from __main__ import fast_gnp_erdos, N, p
import numpy as np
import graph_tool.all as gt
from fibration_mains import FFPartitioning, MBColoring

g = fast_gnp_erdos(N, p, gdirected=True)
fiber_index = g.new_vertex_property('int')
regulation = g.new_edge_property('int')
for n in g.edges(): regulation[n] = 0
g.edge_properties['regulation'] = regulation
g.vertex_properties['fiber_index'] = fiber_index
 ''' 

    TEST_CODE = '''
FFPartitioning(g) '''

    times = timeit.repeat(setup=SETUP_CODE, 
                          stmt=TEST_CODE, 
                          repeat=2, 
                          number=n_repeat)

    print('{}'.format(min(times)/n_repeat))
############################################

####### MINIMAL BALANCED COLORING ########
def MBC_time():
    SETUP_CODE = '''
from __main__ import fast_gnp_erdos, N, p
import numpy as np
import graph_tool.all as gt
from fibration_mains import FFPartitioning, MBColoring

g = fast_gnp_erdos(N, p, gdirected=True)
regulation = g.new_edge_property('int')
for n in g.edges(): regulation[n] = 0
g.edge_properties['regulation'] = regulation'''

    TEST_CODE = '''
MBColoring(g)   '''

    times = timeit.repeat(setup=SETUP_CODE, 
                          stmt=TEST_CODE, 
                          repeat=2, 
                          number=n_repeat)

    print('{}'.format(min(times)/n_repeat))
############################################

if __name__ == "__main__":
    if mode == "ffp":   FFP_time()
    elif mode == "mbc": MBC_time()
    else: print("flags: 'fpp', 'mbc'")