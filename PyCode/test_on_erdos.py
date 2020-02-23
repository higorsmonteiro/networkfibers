import sys
from utils import *
import graph_tool.all as gt
import matplotlib.pyplot as plt

N = [512, 1024, 2048, 4096, 8192]
k_aver = [0.5, 1., 2., 4., 8., 16., 32.]

ps = []
giant_size = []
dir_giant_size = []

for index, n in enumerate(N):
	ps.append([])
	giant_size.append([])
	dir_giant_size.append([])
	for k in k_aver:
		p = k/(n-1)
		ps[index].append(p)
		random_g = fast_gnp_erdos(n, p, gdirected=True)
		dir_giant = gt.extract_largest_component(random_g, directed=True)
		giant = gt.extract_largest_component(random_g, directed=False)

		giant_size[index].append(giant.num_vertices()/n)
		dir_giant_size[index].append(dir_giant.num_vertices()/n)

num_row = 1
num_col = 2
fig, (AX1, AX2) = plt.subplots(num_row, num_col, figsize=(10,6))

for index, n in enumerate(N):
	AX1.plot(k_aver, giant_size[index], marker="o", label=str(n))
	AX2.plot(k_aver, dir_giant_size[index], marker="s", label=str(n))

plt.show()

## Conclusion, the network is pretty much connected with k_aver larger 
## than 4(weak and strong components).







