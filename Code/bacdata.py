import numpy as np 
import pandas as pd
from collections import defaultdict


name_relation = defaultdict(list)
#table = pd.read_csv("../Data/rawData/interactions_bacilus.txt", sep=" ", names=['Source', 'Target'])
table = pd.read_csv("../Data/Ecoli/Ecoli.txt", sep=" ", names=['Source', 'Target', 'Regulation'])

destname = "../Data/ECOLINgenes.dat"
edgelistname = "../Data/ECOLIedgelist.dat"

number_lines = len(table.index)
edgelist = np.zeros([number_lines, 2], int)

sources = list(table['Source'][:])
targets = list(table['Target'][:])
types = list(table['Regulation'][:])
sources = sources + targets
sources = list(set(sources))

print(set(types))

ngenes = len(sources)
labels = np.arange(0, ngenes, 1)

type_reg = []
for index, row in table.iterrows():
    #if index==0: continue
    src = row['Source']
    dest = row['Target']
    edgelist[index, 0] = labels[sources.index(src)]
    edgelist[index, 1] = labels[sources.index(dest)]
    type_reg.append(types[index])

file = open(destname, "w")
file.write("%d\n" % ngenes)
file.close()

file = open(edgelistname, "w")
for index, type_str in enumerate(type_reg):
    file.write("%d\t%d\t%s\n" % (edgelist[index,0], edgelist[index,1], type_str))
file.close()