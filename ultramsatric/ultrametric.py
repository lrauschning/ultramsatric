from typing import Set
import tempfile

import dendropy
import numpy as np

from .distance import DistMat

class Tree:
    def __init__(self, l, ldist: float, r, rdist: float):
        self.l = l
        self.ldist = ldist
        self.r = r
        self.rdist = rdist

    def __init__(self, l, r):
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

def UPGMA(d: DistMat) -> Tree:
    """
    Implementation of the UPGMA algorithm
    or use from scipy/dendropy?
    """
    #TODO same with NJ/minimal evolution
    #TODO create class for tree/use dendropy? or keep as string in newick format?
    pass

def UPGMA_matrix(d: DistMat) -> DistMat:
    """
    Quick & Dirty function to get a UPGMA distance matrix from a distance matrix by calling the UPGMA implementation in dendropy directly.
    Does not use the Tree class.
    """
    tmp = tempfile.NamedTemporaryFile(mode='wt').name + '.tsv'
    #print(tmp)
    d.to_dendropy_csv(tmp)
    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(tmp, mode='rt'), delimiter='\t')
    return DistMat.from_dendropy(pdm.upgma_tree().phylogenetic_distance_matrix())



def closest_ultrametric(matrix: DistMat) -> str:
    pass

def minimal_spanning_tree(matrix: DistMat) -> str:
    # implementable in linear time for a dense distance
    pass

