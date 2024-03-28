#!/bin/python3

from typing import List, Callable, Dict
import os
import itertools
import math

import numpy as np
import dendropy

from .msa import MSA
from .substitutions import *


def alndist(ref: List[chr], alt: List[chr], subs: Callable[[chr, chr], float] = identity, gapcost: Callable[[int], float] = linear, match_gaps=False) -> float:
    """This function calculates the alignment distance between `ref` and `alt`, using the specified substitution model `subs` and the gapcost function `gapcost`.
    `match_gaps` specifies whether to calculate substitution costs between gaps and residues; in that case, `subs` must accept being called with '-' as a gap symbol and a residue.
    :returns: Alignment distance.
    """
    if not len(ref) == len(alt):
        raise ValueError("The sequences must have equal length!")

    refit = iter(ref)
    altit = iter(alt)
    gaplen = 0
    dist = 0
    try:
        refch = next(refit)
        altch = next(altit)
        while refit or altit:
            ## check for insertion
            while refch == '-':
                if altch == '-': # ignore gaps in both sequences
                    #if match_gaps: # uncomment to allow gap-gap penalties >0
                    #    dist += subs('-', '-')
                    refch = next(refit)
                    altch = next(altit)
                    continue
                if match_gaps: # align residue to gap if specified
                    dist += subs(refch, altch)
                gaplen += 1
                refch = next(refit)
            if gaplen > 0:
                dist += gapcost(gaplen)
                gaplen = 0
                continue

            ## check for deletion
            while altch == '-':
                if refch == '-': # ignore gaps in both sequences
                    #if match_gaps: # uncomment to allow gap-gap penalties >0
                    #    dist += subs('-', '-')
                    refch = next(refit)
                    altch = next(altit)
                    continue
                if match_gaps: # align residue to gap if specified
                    dist += subs(refch, altch)
                gaplen += 1
                altch = next(altit)
            if gaplen > 0:
                dist += gapcost(gaplen)
                gaplen = 0
                continue

            ## must be match
            #print(refch, altch, dist, subs(refch, altch))
            dist += subs(refch, altch)

            ## compare the next positions
            refch = next(refit)
            altch = next(altit)
    except StopIteration:
        # one or both sequences have ended; treat the rest as a gap
        reflen = len([x for x in refit])
        altlen = len([x for x in altit])
        if reflen + altlen + gaplen > 0:
            dist += gapcost(reflen + altlen + gaplen)

    return dist

def log_alndist(ref, alt, subs: Callable[[chr, chr], float] = blosum, gapcost: Callable[[int], float] = affine) -> float:
    return math.log(alndist(ref, alt, subs=subs, gapcost=gapcost))

def sq_alndist(ref, alt, subs: Callable[[chr, chr], float] = blosum, gapcost: Callable[[int], float] = affine) -> float:
    return alndist(ref, alt, subs=subs, gapcost=gapcost)**2

def scoredist(ref:List[chr], alt:List[chr], gapcost=no_gaps, subs=blosum) -> float:
    """Implements the Scoredist protein distance function, as described in https://doi.org/10.1186/1471-2105-6-108.
    Uses the BLOSUM62 matrix and no gap penalty by default.
    Gapcost and substitution costs can be configured using the `gapcost` and `subs` parameters.
    :returns: Scoredist distance between `ref` and `alt`.
    """
    c = 1.3370 # from the paper

    # get expected value of substitution matrix
    #aas = [x for x in BLOSUM][:-5] # remove values coding for unknown AAs
    #ev = sum(map(lambda x: x[0][1]*x[1][1]*BLOSUM[x[0][0]][x[1][0]],
    #    itertools.product(AA_FREQS.items(), AA_FREQS.items()))) / (len(aas)**2)

    #ev = -23.42 # the code block above computes to this

    ev = get_ev(subs)
    #print(ev)

    l = max(len(ref), len(alt)) # get alignment length

    dist = alndist(ref, alt, subs=subs, gapcost=gapcost)
    normdist = max(1, dist - l*ev) # normalize by sequence length, enforce distance >= 1 as a pseudocount in case of above-random dissimilarity
    
    lim = (alndist(ref, ref, subs=subs, gapcost=None) +
           alndist(alt, alt, subs=subs, gapcost=None)) / 2 # there shouldn't be any gaps aligning a sequence to itself, so do not set gapcost
    normlim = max(1, lim - l*ev) # pseudocount again
    #print(l, ev, dist, normdist, lim, normlim)
    
    return -c*math.log(normdist / normlim)*100 # logtransform, scale and return

class DistMat:
    """Class representing a distance matrix.
    The underlying representation is a linearization of an upper triangle matrix lacking the diagonal (as it will always be 0).
    A matrix M = [
    [-, 0, 1, 2, 3],
    [-, -, 4, 5, 6],
    [-, -, -, 7, 8],
    [-, -, -, -, 9],
    [-, -, -, -, -]
    ]
    will be linearized to [0, 1, 2, 3, 4, 5, 6, 7, 8, 9].
    The associated _index method will return the 1-dimensional index in the linearized representation corresponding to the 2-dimensional index in the distance matrix.
    The dimension of the distance matrix is stored at `self.n`.
    """
    def __init__(self, n: int, idmap: Dict[str, int], backing: np.ndarray):
        self.n = n # number of sequences stored
        self.idmap = idmap # map storing the index of each FASTA ID
        self._backing = backing # an upper triangle matrix lacking the diagonal, linearized to a 1D-Array

    @classmethod
    def from_dendropy(cls, pdm: dendropy.PhylogeneticDistanceMatrix):
        n = len(pdm.taxon_namespace)
        idmap = {str(t):id for id, t in enumerate(sorted(pdm.taxon_namespace))}
        backing = np.ndarray(n*(n-1)//2, dtype=np.float32)
        for i in range(n):
            for j in range(i+1, n):
                t = pdm.taxon_namespace[i]
                v = pdm.taxon_namespace[j]
                backing[DistMat.index(idmap[str(t)], idmap[str(v)], n)] = pdm.distance(t, v)

        return cls(n, idmap, backing)
        
    def _index(self, a:int, b:int) -> int:
        return DistMat.index(a, b, self.n)

    def _revindex(self, x:int) -> int:
        return DistMat.revindex(x, self.n)

    def get(self, a:str, b:str) -> float:
        return self._get(self.idmap[a], self.idmap[b])

    def _get(self, a: int, b: int) -> float:
        if b < a: # ensure a <= b
            a, b = b, a
        elif a == b:
            return 0
        return self._backing[self._index(a, b)]

    def set(self, a:str, b:str, v:float):
        self._set(self.idmap[a], self.idmap[b])

    def _set(self, a:int, b:int, v:float):
        self._backing[self._index(a, b)] = v

    def __len__(self) -> int:
        return self.n


    def apply(self, fun):
        """Takes a function taking as arguments the position and current value, and stores the result of applying that function at each position in the Distance Matrix.
        """
        for i in range(self.n):
            for j in range(i+1, self.n):
                #print("old:", self._backing[self._index(i, j)])
                #print("new:", fun(i, j, self._get(i, j)))
                self._backing[self._index(i, j)] = fun(i, j, self._get(i, j))


    @classmethod
    def index(cls, a:int, b:int, n:int) -> int:
        if b < a: # ensure a <= b
            a, b = b, a
        return (a * (a-1))//2 + a * (n - a) + (b - a - 1)

    @classmethod
    def revindex(cls, x:int, n:int) -> (int, int):
        """
        Inverse of the index function. This is unique because indices are discrete.
        index(*revindex(i, n), n) will hold for any i and n.
        """
        a: int = int((2*n-1-math.sqrt(4*n*(n-1)-8*x+1))//2) # look for maximal a, then shift leftover into b
        b: int = int(x + 1 - a*(2*n - 3 - a)/2)
        assert(a <= b)
        return a, b

    def to_full_matrix(self, rnd: int=-1) -> np.ndarray:
        """
        Returns a full n*n matrix containing the distances encoded in this matrix.
        The returned matrix will be symmetric.
        :args: rnd=0: Number of digits to round the returned matrix to, default -1 (no rounding).
        Mainly intended for displaying output graphically.
        """
        ret = np.zeros([self.n, self.n], dtype=self._backing.dtype)
        for i in range(self.n):
            for j in range(self.n):
                ret[(i, j)] = self._get(i, j) if rnd < 0 else round(self._get(i, j), rnd)
        return ret

    def to_dendropy_csv(self, path: os.PathLike, sep='\t', newline='\n'):
        ids = sorted(self.idmap.keys())
        with open(path, 'wt') as f:
            f.write(sep.join([''] + ids))
            f.write(newline)
            for i in range(len(ids)):
                f.write(sep.join([ids[i]] + [str(round(self._get(i, j), 2)) for j in range(len(ids))]))
                f.write(newline)


    def __repr__(self) -> str:
        return "\n".join(["\t".join(map(str, x[:])) for x in self.to_full_matrix(rnd=2)[:]])

    @classmethod
    def from_msa(cls, m: MSA, distfun):
        # init variables
        ids = sorted(m.alns.keys())
        n = len(ids)
        #print(n, ids)
        # stolen from https://stackoverflow.com/a/1679702
        idmap = dict(map(reversed, enumerate(ids)))
        backing = np.ndarray(n*(n-1)//2, dtype=np.float32)

        # calculate pairwise distances
        for i in range(len(ids)):
            for j in range(i+1, len(ids)):
                #print("comparing", ids[i], ids[j])
                backing[DistMat.index(i, j, n)] = distfun(m.alns[ids[i]], m.alns[ids[j]])

        return cls(n, idmap, backing)

    def __sub__ (self, other):
        if self.n != other.n:
            raise ValueError("Tried to subtract Matrices with different dimensions!")
        #if self.idmap != other.idmap:
        #    raise ValueError("Tried to subtract Matrices with different taxons!")
        return DistMat(self.n, self.idmap, self._backing - other._backing)

    ## Implement a few Matrix norms
    def abssum(self) -> float:
        return np.sum(np.abs(self._backing))

    def frobenius(self) -> float:
        return math.sqrt(np.sum(self._backing**2))

    def norm_frobenius(self) -> float:
        return self.frobenius()/len(self._backing)

    def absavg(self) -> float:
        return self.abssum()/len(self._backing)

    def corr(self, other) -> float:
        assert(len(self) == len(other))
        cov = np.sum((self._backing - np.mean(self._backing))*(other._backing - np.mean(other._backing)))/(self.n**2)
        return cov/(np.std(self._backing)*np.std(other._backing))

if __name__ == "__main__":
    import sys
    ref = [x for x in sys.argv[1]]
    alt = [x for x in sys.argv[2]]
    print(alndist(ref, alt))
