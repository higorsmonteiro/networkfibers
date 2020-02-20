import numpy as np
from utils import *
from fiber import *
import FFPf as ffp
import graph_tool.all as gt
import minimalcoloringf as mbc
from collections import Counter, deque

def MBColoring(g):
    N = g.get_vertices()

    # Properties of the nodes: colors and ISCV.
    node_colors = g.new_vertex_property('int')
    iscv = g.new_vertex_property('string')
    g.vertex_properties['node_colors'] = node_colors
    g.vertex_properties['iscv'] = iscv

    #### INITIALIZATION: Criterion -> inputless SCC's as different classes. ####
    fibers = mbc.Initialization(g)  # List of fiber classes.
    mbc.set_colors(g, fibers)       # Set the colors for each node according its fiber.

    ncolor_after = len(fibers)
    ncolor_before = -1

    mbc.set_ISCV(g, ncolor_after)

    ######### REFINEMENT LOOP ############
    while ncolor_after!=ncolor_before:
        iscv_list = list(g.vp.iscv)
        splitted = []
        ''' For each fiber, we split it according the value of
            the ISCV of each node inside the fiber. '''
        for class_index, fblock in enumerate(fibers):
            if fblock.get_number_nodes() <= 1: continue

            # defines the list of nodes of the current fiber and their ISCVs.
            fiber_nodeindex = [node for node in fblock.get_nodes()]
            fiber_iscv = [iscv_list[node] for node in fblock.get_nodes()]
            iscv_count = Counter(fiber_iscv)
            if len(iscv_count) == 1: continue # The fiber is not splitted.

            # Now we split the fiber according their ISCV 'fiber_iscv'.
            splitted_list = mbc.split_fiberf(class_index, fblock, fiber_nodeindex, fiber_iscv)

            # If the current fiber is splitted, remove the father from the list.
            #fibers.remove(fblock)
            for new_fiber in splitted_list: fibers.append(new_fiber)

        ncolor_before = ncolor_after
        ncolor_after = len(fibers)
        mbc.set_colors(g, fibers)

        mbc.set_ISCV(g, ncolor_after)


def FFPartitioning(g):
    bqueue = deque([])
    partition = ffp.Initialization(g, bqueue)

    # Until the queue is empty, we procedure the splitting process.
    while bqueue:
        pivot_set = bqueue.popleft()
        ffp.input_splitf(partition, pivot_set, g, 1, bqueue)

if __name__=="__main__":
    N = 32 
    p = 1/(N-1)
    
    g = fast_gnp_erdos(N, p, gdirected=True)
    print(g.num_vertices(), g.num_edges())
    regulation = g.new_edge_property('int')
    for n in g.edges(): regulation[n] = 0
    g.edge_properties['regulation'] = regulation

    edge_color = g.new_edge_property('int')
    color_name = g.new_vertex_property('string')
    g.edge_properties['edge_color'] = edge_color
    g.vertex_properties['color_name'] = color_name
    
    MBColoring(g)
    node_colors = g.vp.node_colors
    for v in g.get_vertices(): color_name[v] = str(node_colors[v])
    for e in g.edges():
        source = e.source()
        edge_color[e] = node_colors[source]

    node_colors = g.vp.node_colors
    gt.graph_draw(g, vertex_text=color_name, vertex_color=node_colors, edge_color=edge_color, vertex_fill_color=node_colors, output_size=(800, 800), output="../Figures/fibers/test_mbc.pdf")
