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