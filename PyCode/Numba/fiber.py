import numpy as np
from numba import jitclass, int64, float64

spec = [
    ('index', int64),
    ('number_nodes', int64),
    ('branching', float64),
    ('reg_number', int64),
    ('regulators', int64),
    ('fibernodes', int64)
]

@jitclass(spec)
class FiberBlock:
    def __init__(self):
        self.index = -1
        #self.regtype = [-1, -1, -1]
        self.number_nodes = 0
        self.branching = -1.0
        self.reg_number = -1
        self.regulators = [int64(x) for x in range(0)]
        self.fibernodes = [int64(x) for x in range(0)]

    def insert_node(self, node):
        self.number_nodes += 1
        self.fibernodes.append(node)
    
    def delete_nodes(self, nodelist):
        self.number_nodes -= len(nodelist)
        self.fibernodes = [item for item in self.fibernodes if item not in nodelist]
    
    #def delete_node(self, node):
    #    self.number_nodes -= 1
    #    self.fibernodes.remove(node)

    def insert_nodelist(self, nodelist):
        self.number_nodes += len(nodelist)
        self.fibernodes += nodelist

    #def insert_regulator(self, reg):
    #    self.number_reg += 1
    #    self.regulators.append(reg)

    def get_nodes(self):
        return self.fibernodes

    def get_number_nodes(self):
        return self.number_nodes

    def get_regulators(self):
        return self.regulators

    def show_nodes(self):
        print(self.fibernodes)

    def sucessor_nodes(self, graph):
        ''' Lists all nodes that receives information
            from 'self' fiber. This listing is useful
            for pivot sets, where we need to find the
            unstable classes with respect to 'self'. '''
        sucessors = []
        for node in self.fibernodes:
            out_neigh = graph.get_out_neighbors(node)
            sucessors += list(out_neigh)

        return list(set(sucessors))


    def show_nodes_name(self, graph):
        '''
            Used by the 'PrintFibers' function.
        '''
        names = graph.vertex_properties['node_names']
        for node in self.fibernodes:
            print(names[node], end=" ")
        print("")

    def input_stability(self, graph, Set, regulation):
        '''
            Given a fiber 'Set' and the graph with its regulation
            types, then it checks if the fiber is input-set stable
            with respect to 'Set'.
        '''
        set_nodes = Set.get_nodes()
        edges_received = np.zeros([3, len(self.fibernodes)], int)
        
        for setnode in set_nodes:
            out_edges = graph.get_out_edges(setnode, [graph.edge_index])
            for edge in out_edges:
                if edge[1] in self.fibernodes:
                    cur_index = self.fibernodes.index(edge[1])
                else: continue
                edges_received[regulation[edge[2]], cur_index] += 1

        col = edges_received[:,0]
        for k in range(len(self.fibernodes)):
            if not np.array_equal(col, edges_received[:,k]):
                return -1
        return 1


class StrongComponent:
    def __init__(self):
        self.nodes = []
        self.have_input = False
        self.number_nodes = 0
        self.type = -1

    def insert_node(self, node):
        self.number_nodes += 1
        self.nodes.append(node)

    def get_nodes(self):
        return self.nodes

    def show_nodes(self):
        print(self.nodes)

    def check_input(self, graph):
        '''
            Check if the SCC receives or not 
            input from other components of the
            network.
        '''
        for v in self.nodes:
            input_neigh = graph.get_in_neighbors(v)
            for neigh in input_neigh:
                try: bl = self.nodes.index(neigh)
                except: bl = -1
                # the SCC have input from other component.
                if bl == -1:    
                    self.have_input = True
                    break
            if self.have_input==True: break

    def classify_strong(self, graph):
        ''' It should be called after 'check_input' '''
        
        # if the SCC has input, then it belongs to the main fiber.
        if self.have_input: self.type = 0
        else:
            #   If it doesn't receive any external input, then we
            #   must check if it is an isolated autorregulated
            #   node.
            if self.number_nodes==1:
                input_neigh = graph.get_in_neighbors(self.nodes[0])
                if input_neigh.shape[0] == 0:
                    self.type = 1
                else: self.type = 2 # Isolated autorregulated node. 
            else: self.type = 1     # SCC does not receive any external input.


    def show_input_bool(self):
        print(self.have_input)

    def get_input_bool(self):
        return self.have_input 
