import numpy as np

class FiberBlock:
    def __init__(self):
        self.index = -1
        self.regtype = [-1, -1, -1]
        self.number_nodes = 0
        self.branching = -1.0
        self.reg_number = -1
        self.regulators = []
        self.fibernodes = []

    def insert_node(self, node):
        self.number_nodes += 1
        self.fibernodes.append(node)

    def insert_regulator(self, reg):
        self.number_reg += 1
        self.regulators.append(reg)

    def get_nodes(self):
        return self.fibernodes

    def get_number_nodes(self):
        return self.number_nodes

    def get_regulators(self):
        return self.regulators

    def show_nodes(self):
        print(self.fibernodes)

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
