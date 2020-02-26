import sys
from utils import *
import graph_tool.all as gt
from main_f import MBColoring
from collections import Counter
import matplotlib.pyplot as plt
from minimalcoloringf import get_possible_unstable_classes, fast_checking

N = int(sys.argv[1])
k_aver = float(sys.argv[2])
p = k_aver/(N-1)

random_g = fast_gnp_erdos(N, p, gdirected=True)
fibers = MBColoring(random_g, get_flist=True)

print("# of fibers = %i; N = %i; p = %lf; <k> = %lf" %(len(fibers), N, p, k_aver))

n_edgetype = 1
for pivot in fibers:
	eta = pivot.sucessor_nodes(random_g)
	regulation = random_g.edge_properties['regulation'].a
    # Given the node number, 'f' gives its index in 'eta'.
	f = defaultdict(lambda:-1)
	for eta_index, sucessor in enumerate(eta):  f[sucessor] = eta_index

    # 'R' represents a matrix (n_edgetype, len(eta)).
	R = np.vstack([np.zeros(len(eta), int) for row in range(n_edgetype)])

	calc_R(R, random_g, pivot, f, regulation)
	receiver_classes = get_possible_unstable_classes(random_g, pivot, fibers)
	fast_checking(receiver_classes, eta, f, R, fibers, n_edgetype, random_g)

#PrintFibers(fibers, random_g)
