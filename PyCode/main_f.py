import sys
import numpy as np
import FFPf as ffp
import MBCf as mbc
from utils import *
from fiber import *
import graph_tool.all as gt
from collections import Counter, deque

def MBColoring(g, get_flist=False, num_edgetype=1):
    '''
        Given a network 'g', this function is responsible
        to partition the set of nodes into classes where nodes
        are input-tree symmetrical.

        The network 'g' must have an edge property called 
        'regulation' having integer values corresponding to
        the type of all edges.

        if 'get_flist' is set TRUE, then this function returns
        the list of fibers resulted from the partitioning.
    '''

    ########### Properties of the network ############
    fiber_index = g.new_vertex_property('int')  # colors
    iscv = g.new_vertex_property('string')
    edge_color = g.new_edge_property('int')   
    color_name = g.new_vertex_property('string')    # color index string.
    regulation = g.new_edge_property('int')         # type of edge.
    g.vertex_properties['fiber_index'] = fiber_index
    g.vertex_properties['iscv'] = iscv
    g.edge_properties['edge_color'] = edge_color
    g.vertex_properties['color_name'] = color_name
    g.edge_properties['regulation'] = regulation
    ##################################################

    #### INITIALIZATION: Criterion -> inputless SCC's as different classes. ####
    fibers = mbc.Initialization(g)  # List of fiber classes.
    ncolor_after = len(fibers)
    ncolor_before = -1
    # Set the colors for each node according its 'fibers' index.
    mbc.set_colors(g, fibers)
    ############################################################################       

    ######### REFINEMENT LOOP ############
    mbc.set_ISCV(g, ncolor_after, num_edgetype)
    while ncolor_after!=ncolor_before:
        iscv_list = list(g.vp.iscv)
        ''' For each fiber, we split it according the value of
            the ISCV of each node inside the fiber. '''
        for class_index, fclass in enumerate(fibers):
            if fclass.get_number_nodes() <= 1: continue
            # defines the list of nodes of the current fiber and their ISCVs.
            f_nodeindex = [node for node in fclass.get_nodes()]
            fiber_iscv = [iscv_list[node] for node in fclass.get_nodes()]
            iscv_count = Counter(fiber_iscv)
            if len(iscv_count) == 1: continue # The fiber is not splitted.
            # Now we split the fiber according their ISCV 'fiber_iscv'.
            splitted_list = mbc.split_fiberf(class_index, fclass, f_nodeindex, fiber_iscv)
            for new_fiber in splitted_list: fibers.append(new_fiber)
        ncolor_before = ncolor_after
        ncolor_after = len(fibers)
        mbc.set_colors(g, fibers)

        mbc.set_ISCV(g, ncolor_after, num_edgetype)

    # Properties for network drawing.
    for v in g.get_vertices(): color_name[v] = str(fiber_index[v])
    for e in g.edges():
        source = e.source()
        edge_color[e] = fiber_index[source]

    if get_flist==True: return fibers
    

def FFPartitioning(g, get_flist=False, num_edgetype=1):
    # Necessary property maps.
    edge_color = g.new_edge_property('int')
    color_name = g.new_vertex_property('string')
    fiber_index = g.new_vertex_property('int')
    regulation = g.new_edge_property('int')
    g.edge_properties['edge_color'] = edge_color
    g.vertex_properties['color_name'] = color_name
    g.vertex_properties['fiber_index'] = fiber_index
    g.edge_properties['regulation'] = regulation
    ###################################################
    for n in g.edges(): regulation[n] = 0   # One edge type.

    bqueue = deque([])
    partition = ffp.Initialization(g, bqueue)

    # Until the queue is empty, we procedure the splitting process.
    while bqueue:
        pivot_set = bqueue.popleft()
        ffp.input_splitf(partition, pivot_set, g, num_edgetype, bqueue)

    if get_flist==True: return partition

if __name__=="__main__":
    N = int(sys.argv[1])
    k_aver = float(sys.argv[2])
    p = k_aver/(N-1)
    
    nedgetype = 200
    g = fast_gnp_erdos(N, p, num_edgetype=nedgetype, gdirected=True)
    f = MBColoring(g, num_edgetype=nedgetype, get_flist=True)
    print(len(f))
    #draw_fibers(g, output_filename="test.pdf", edgetext=True, layout='spring')

    
    
