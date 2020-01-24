import numpy as np
import graph_tool.all as gt
from collections import Counter

def get_SOLITAIRE(graph):
    solitaires = []
    for v in graph.get_vertices():
        if (graph.get_in_neighbors(v)).shape[0] == 0:
            solitaires.append(v)
    return solitaires

def set_SOLITAIRE(graph):
    solitaires = []
    for v in graph.get_vertices():
        if (graph.get_in_neighbors(v)).shape[0] == 0:
            solitaires.append(v)

    node_colors = graph.vp.node_colors
    for node in graph.get_vertices():
        node_colors[node] = 0
    
    label = 1
    for sol in solitaires:
        node_colors[sol] = label
        label += 1
    return len(solitaires) + 1

# Define the ISCV for each node. At this step, each ISCV has size 'ncolor'.
def set_ISCV(graph, ncolor):
    # intrinsic properties
    iscv = graph.vp.iscv
    node_colors = graph.vp.node_colors

    # initializate iscv
    for v in graph.get_vertices():
        iscv[v] = ""

    # define the iscv for each node.
    for v in graph.get_vertices():
        in_neighbors = graph.get_in_neighbors(v)
        input_colors = []
        for neigh in in_neighbors:
            input_colors.append(node_colors[neigh])
        # Counts how many inputs the 'node' receives from each color.
        Colors_counter = Counter(input_colors)
        for k in range(ncolor):
            # 'Color_counter[k] returns the number of inputs from color 'k'.
            iscv[v] += str(Colors_counter[k])


def print_colors(graph):
    node_colors = graph.vp.node_colors
    for v in graph.get_vertices():
        print(v, node_colors[v])

def number_colors(graph):
    node_colors = list(graph.vp.node_colors.a)
    print(len(set(node_colors)))
