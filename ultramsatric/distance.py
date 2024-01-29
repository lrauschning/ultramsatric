#!/bin/python3

from typing import List, Callable, Dict
import os
import itertools
import math

import numpy as np
import dendropy
import blosum as bl

from .msa import MSA


BLOSUM = bl.BLOSUM(62)

AA_FREQS = {'A': 8.76, 'R': 5.78, 'N': 3.93, 'D': 5.49, 'C': 1.38, 'Q': 3.9,
           'E': 6.32, 'G': 7.03, 'H': 2.26, 'I': 5.49, 'L': 9.68, 'K': 5.19, 'M': 2.32,
           'F': 3.87, 'P': 5.02, 'S': 7.14, 'T': 5.53, 'W': 1.25, 'Y': 2.91, 'V': 6.73
        }
# from https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties


def identity(ref:chr, alt:chr) -> float:
    return 1 if ref != alt else 0

def blosum(ref:chr, alt:chr) -> float:
    return BLOSUM[ref][alt]

def linear(n:int) -> float:
    return 3*n

def affine(n:int) -> float:
    return 3 + 2*n

def no_gaps(n:int) -> float:
    return 0.0

def alignment_distance(ref: List[chr], alt: List[chr], subs: Callable[[chr, chr], float] = identity, gapcost: Callable[[int], float] = linear) -> float:
    """This function calculates the alignment distance between `ref` and `alt`, using the specified substitution model `subs` and the gapcost function `gapcost`.
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
                    refch = next(refit)
                    altch = next(altit)
                    continue
                gaplen += 1
                refch = next(refit)
            if gaplen > 0:
                dist += gapcost(gaplen)
                gaplen = 0
                continue

            ## check for deletion
            while altch == '-':
                if refch == '-': # ignore gaps in both sequences
                    refch = next(refit)
                    altch = next(altit)
                    continue
                gaplen += 1
                altch = next(altit)
            if gaplen > 0:
                dist += gapcost(gaplen)
                gaplen = 0
                continue

            ## must be match
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

def scoredist(ref:List[chr], alt:List[chr], gaps=False, blo=62) -> float:
    """Implements the Scoredist protein distance function, as described in https://doi.org/10.1186/1471-2105-6-108.
    Uses the BLOSUM62 matrix and no gap penalty by default.
    The BLOSUM matrix can be set as the `blo` parameter.
    If `gaps` is set to `True`, it uses an affine gap penalty (not recommended in the paper).
    :returns: Scoredist distance between `ref` and `alt`.
    """
    c = 1.3370 # from the paper

    BLOSUM = bl.BLOSUM(blo)

    # get expected value of substitution matrix
    #aas = [x for x in BLOSUM][:-5] # remove values coding for unknown AAs
    #ev = sum(map(lambda x: x[0][1]*x[1][1]*BLOSUM[x[0][0]][x[1][0]],
    #    itertools.product(AA_FREQS.items(), AA_FREQS.items()))) / (len(aas)**2)

    ev = -23.42 # the code block above computes to this

    l = max(len(ref), len(alt)) # get alignment length

    alndist = alignment_distance(ref, alt, subs=blosum, gapcost= (lambda n: -8 -3*n) if gaps else (lambda n: -4*n))
    normdist = max(1, alndist - l*ev) # normalize by sequence length, enforce distance >= 1 as a pseudocount in case of above-random dissimilarity
    
    lim = (alignment_distance(ref, ref, subs=blosum, gapcost=None) +
           alignment_distance(alt, alt, subs=blosum, gapcost=None)) / 2 # there shouldn't be any gaps aligning a sequence to itself, so do not set gapcost
    normlim = max(1, lim - l*ev) # pseudocount again
    #print(l, ev, alndist, normdist, lim, normlim)
    
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
        

    def get(self, a:str, b:str) -> float:
        return self._get(self.idmap[a], self.idmap[b])

    def _get(self, a: int, b: int) -> float:
        if b < a: # ensure a <= b
            a, b = b, a
        elif a == b:
            return 0
        return self._backing[self._index(a, b)]

    def _index(self, a:int, b:int) -> int:
        return DistMat.index(a, b, self.n)

    def apply(self, fun):
        """Takes a function taking as arguments the position and current value, and stores the result of applying that function at each position in the Distance Matrix.
        """
        for i in range(len(ids)):
            for j in range(i+1, len(ids)):
                self._backing[sevf._index(i, j)] = fun(i, j, self._get(i, j))


    @classmethod
    def index(cls, a:int, b:int, n:int) -> int:
        if b < a: # ensure a <= b
            a, b = b, a
        return (a * (a-1))//2 + a * (n - a) + (b - a - 1)

    def to_full_matrix(self) -> np.ndarray:
        ret = np.zeros([self.n, self.n], dtype=self._backing.dtype)
        for i in range(self.n):
            for j in range(self.n):
                ret[(i, j)] = self._get(i, j)
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
        return str(self.to_full_matrix())

    @classmethod
    def from_msa(cls, m: MSA, distfun):
        # init variables
        ids = sorted(m.alns.keys())
        n = len(ids)
        print(n, ids)
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

if __name__ == "__main__":
    import sys
    ref = [x for x in sys.argv[1]]
    alt = [x for x in sys.argv[2]]
    print(alignment_distance(ref, alt))
