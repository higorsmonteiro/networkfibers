'''
    Based on the methods presented on the draft, we collect the
    runtime measurements for each fibration algorithm and we store
    on data files.
'''
import sys
import timeit
import numpy as np
import graph_tool.all as gt

N = int(sys.argv[1])
k = float(sys.argv[2])
mode = sys.argv[3]
n_edgetype = int(sys.argv[4])

p = k/(N-1)
n_repeat = 10

######## FAST FIBRATION PARTITIONING ########
def FFP_time():
    SETUP_CODE = '''
from __main__ import N, p, n_edgetype
from utils import fast_gnp_erdos
import numpy as np
import graph_tool.all as gt
from main_f import FFPartitioning

g = fast_gnp_erdos(N, p, num_edgetype=n_edgetype, gdirected=True)
''' 

    TEST_CODE = '''
FFPartitioning(g, num_edgetype=n_edgetype) '''

    times = timeit.repeat(setup=SETUP_CODE, 
                          stmt=TEST_CODE, 
                          repeat=1, 
                          number=n_repeat)

    print(times[0]/n_repeat)
############################################

####### MINIMAL BALANCED COLORING ########
def MBC_time():
    SETUP_CODE = '''
from __main__ import N, p, n_edgetype
from utils import fast_gnp_erdos, in_degree_average, out_degree_average
import numpy as np
import graph_tool.all as gt
from main_f import MBColoring

g = fast_gnp_erdos(N, p, num_edgetype=n_edgetype, gdirected=True)
'''

    TEST_CODE = '''
MBColoring(g, num_edgetype=n_edgetype)   '''

    times = timeit.repeat(setup=SETUP_CODE, 
                          stmt=TEST_CODE, 
                          repeat=1, 
                          number=n_repeat)

    print(times[0]/n_repeat)
############################################

if __name__ == "__main__":
    if mode == "ffp":
        FFP_time()
    elif mode == "mbc": 
        MBC_time()
    else: print("flags: 'fpp', 'mbc'")