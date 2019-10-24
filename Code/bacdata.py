import numpy as np 
import pandas as pd
from collections import defaultdict

name_relation = defaultdict(list)
#table = pd.read_csv("../Data/rawData/interactions_bacilus.txt", sep=" ", names=['Source', 'Target'])
table = pd.read_csv("../Data/Ecoli/Ecoli.txt", sep=" ", names=['Source', 'Target', 'Regulation'])

number_lines = len(table.index)
edgelist = np.zeros([number_lines, 2], int)

sources = list(table['Source'][:])
targets = list(table['Target'][:])
sources = sources + targets
sources = list(set(sources))

labels = np.arange(0, len(sources), 1)

for index, row in table.iterrows():
    if index==0: continue
    src = row['Source']
    dest = row['Target']
    edgelist[index, 0] = labels[sources.index(src)]
    edgelist[index, 1] = labels[sources.index(dest)]

print(edgelist.shape)
np.savetxt("../Data/ECOLINgenes.dat", [labels.shape[0]], fmt="%d")
np.savetxt("../Data/ECOLIedgelist.dat", edgelist, fmt="%d")