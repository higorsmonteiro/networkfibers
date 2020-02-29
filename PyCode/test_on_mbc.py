'''
	Minimal balancing coloring code testing. For several average degree values
	and for all possible number of edge types, this code runs the MBC's code and
	after checks if the final classes are all input-set stable.

	arg1 -> size of the random network.
	arg2 -> average degree desired.
'''

import sys
from utils import *
import graph_tool.all as gt
from main_f import MBColoring
from collections import Counter
import matplotlib.pyplot as plt
from MBCf import check_sucessor_stability

N = int(sys.argv[1])
k_aver = float(sys.argv[2])
p = k_aver/(N-1)

nedgetype = 1	# number of edge types.
random_g = fast_gnp_erdos(N, p, num_edgetype=nedgetype, gdirected=True)
fibers = MBColoring(random_g, num_edgetype=nedgetype, get_flist=True)

print("# of fibers = %i; N = %i; p = %lf; <k> = %lf" %(len(fibers), N, p, k_aver))

'''
	For each class in the final partitioning, we check if all
	classes that receives information from 'pivot' are stable.
'''
for pivot in fibers:
	check_sucessor_stability(pivot, fibers, random_g, nedgetype)

