import numpy as np
from distance import DistMat

def UPGMA(matrix: DistMat) -> Tree:
    """
    Implementation of the UPGMA algorithm
    or use from scipy/dendropy?
    """
    #TODO same with NJ/minimal evolution
    #TODO create class for tree/use dendropy? or keep as string in newick format?

def closest_ultrametric(matrix: DistMat) -> str:
    pass

def minimal_spanning_tree(matrix: DistMat) -> str:
    # implementable in linear time for a dense distance
    pass

class Tree:
    def __init__(self, l: Tree, ldist: float, r: Tree, rdist: float):
        self.l = l
        self.ldist = ldist
        self.r = r
        self.rdist = rdist

    def __init__(self, l: Tree, r: Tree):
        self.l = l
        self.ldist = 0.0
        self.r = r
        self.rdist = 0.0

    def __init__(self):
        self.l = None
        self.ldist = 0.0
        self.r = None
        self.rdist = 0.0

    def isleaf(self) -> bool:
        return self.l is None and self.r is None

    def label(self, label: str):
        self.label = label

    def get_domain(self) -> Set[str]:
        """Returns the domain of the node, that is the set of all node labels of beneath this node"""
        ret = set()
        if self.l:
            ret = ret.union(l.get_domain())
        if self.r:
            ret = ret.union(r.get_domain())
        if self.label:
            ret.add(self.label)
        return ret

