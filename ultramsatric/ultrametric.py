from typing import Set, Tuple, List, Dict
import tempfile

import dendropy
import numpy as np
import math

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

def root_ext_add(d: DistMat) -> DistMat:
    """
    Transforms an additive distance matrix into an ultrametric one by rooting along the maximal edge and extending edges incident to leaves.
    This is guaranteed to return an ultrametric matrix if the input is additive.
    If the input matrix is not additive, may return an invalid DistMat.
    """
    amax, bmax = d._revindex(np.argmax(d._backing)) # select root by getting maximal distance
    dmax = d._get(amax, bmax)

    # construct ultrametric matrix
    um = DistMat(d.n, d.idmap, np.ndarray(d.n*(d.n-1)//2, dtype=np.float32))
    um.apply(lambda i, j, v: dmax - (d._get(amax, j) + d._get(amax, i) - d._get(i, j))/2)
    #print(d)
    #print(um)
    return um


def NJ_matrix(d: DistMat) -> DistMat:
    """
    Quick & Dirty function to get a NJ distance matrix from a distance matrix by calling the NJ implementation in dendropy directly.
    Does not use the Tree class.
    """
    tmp = tempfile.NamedTemporaryFile(mode='wt').name + '.tsv'
    #print(tmp)
    d.to_dendropy_csv(tmp)
    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(tmp, mode='rt'), delimiter='\t')
    return DistMat.from_dendropy(pdm.nj_tree().phylogenetic_distance_matrix())


def tallest_ultrametric(d: DistMat) -> DistMat:
    """
    Implement the algorithm for a closest ultrametric tree of a distance matrix from Prof. Volker Heuns lecture script, section 2.7, page 161 (in version 6.28).
    The algorithm is described in Figure 2.66.
    This implementation does not explicitly construct the tree, but directly computes the patristic distances from it at each step.
    The recursion in the algorithm is unfolded into a single loop to facilitate the matrix construction.
    :returns: A DistMat object representing the ultrametric distance matrix corresponding to the tallest ultrametric tree that is compatible to the input distances. These are not required to be additive or ultrametric.
    """
    mst = mst_from_dmat(d)
    um = DistMat(d.n, d.idmap, np.ndarray(d.n*(d.n-1)//2, dtype=np.float32))
    #print(mst)

    parts = [set(range(d.n))] # initialise with the full set of leaves

    while parts: # so long as there are unsplit sets remaining
        part = parts[0]
        #print("Iterating:", part)
        parts = parts[1:]
        assert(len(part) > 1) # this implementation does not explicitly iterate over leaves
        
        # select the edge with the heighest weight from this partition of the MST
        vmax = 0
        imax = 0
        e_wt = 0
        for v in part:
            for i in mst[v]:
                wt = d._get(v, i)
                if wt > e_wt:
                    vmax = v
                    imax = i
                    e_wt = wt
        #print(vmax, imax, e_wt)
        # remove this edge from the mst
        mst[v].remove(i)
        mst[i].remove(v)
        # select the two partitions
        vpart = dft(v, mst)
        ipart = dft(i, mst)
        #print(vpart, ipart)

        # add the new partitions to process splits in them later
        if len(vpart) > 1: # do not add single-item partitions
            parts.append(vpart)
        if len(ipart) > 1: # do not add single-item partitions
            parts.append(ipart)

        # incorporate this split into the matrix
        for v in vpart:
            for i in ipart:
                um._set(v, i, e_wt)

    #done
    return um



def dft(start: int, graph: Dict) -> Set:
    """
    Performs a Depth-First traversal starting from the specified node.
    :returns: The nodes visited in order.
    """
    #TODO benchmark if this is faster than recursion
    stack = [start]
    for x in stack:
        stack.extend(filter(lambda x: x not in stack, iter(graph[x])))
        # this is n^2, think of a more efficient traversal
    return stack[::-1]

def mst_from_dmat(d: DistMat) -> Dict[int, List[int]]:
    """
    Implements the DJP algorithm on a distance matrix.
    Runs in O(n^2) time and O(n) space.
    :returns: A map of adjacency sets corresponding to the MST.
    """
    rem = set(range(d.n)) # set of nodes that are not yet part of the MST

    ## Init: choose the smallet edge
    amin, bmin = d._revindex(np.argmin(d._backing))
    #mst = [(amin, bmin)] # list of edges in the MST
    mst = {amin: {bmin}, bmin:{amin}} # init map of adjacency sets of the MST
    #done = {amin, bmin} # set of nodes already in the MST
    rem -= {amin, bmin}

    while rem: # iterate until no nodes remain
        min_d = math.inf
        out, new = 0, 0
        # iterate through all edges that could be added in this step
        for i in mst.keys():
            for j in rem:
                dist = d._get(i, j)
                # store the minimal one
                if dist < min_d:
                    min_d = dist
                    out, new = i, j

        # add new node to the MST
        mst[out].add(new)
        mst[new] = {out}
        #mst.append((out, new))
        #done.add(new)
        rem.remove(new)

    #print(mst)
    return mst
