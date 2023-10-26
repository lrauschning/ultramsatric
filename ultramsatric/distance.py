#!/bin/python3
from typing import List, Callable, Dict
import numpy as np
from msa import MSA


def identity(ref:chr, alt:chr) -> float:
    return 1 if ref != alt else 0

def blosum(ref:chr, alt:chr) -> float:
    pass
#blosum = #blosum matrix as np.ndarray

def linear(n:int) -> float:
    return 3*n

def affine(n:int) -> float:
    return 3 + 2*n

def alignment_distance(ref: List[chr], alt: List[chr], dfun: Callable[[chr, chr], float] = identity, gapcost: Callable[[int], float] = linear) -> float:
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
            if refch != altch:
                dist += dfun(refch, altch)

            ## compare the next positions
            refch = next(refit)
            altch = next(altit)
    except StopIteration:
        # one or both sequences have ended; treat the rest as a gap
        reflen = len([x for x in refit])
        altlen = len([x for x in altit])
        dist += gapcost(reflen + altlen + gaplen)

    return dist


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

    @classmethod
    def index(cls, a:int, b:int, n:int) -> int:
        if b < a: # ensure a <= b
            a, b = b, a
        return (a * (a-1))//2 + a * (n - a) + (b - a - 1)

    def to_full_matrix(self) -> np.ndarray:
        ret = np.zeros([self.n, self.n], dtype=self._backing.dtype)
        for i in range(self.n):
            for j in range(self.n):
                ret[(i, j)] = self._backing[self._index(i, j)]
        return ret

    def __repr__(self) -> str:
        return str(self.to_full_matrix())

    @classmethod
    def from_msa(cls, m: MSA, distfun):
        # init variables
        ids = sorted(m.alns.keys())
        n = len(ids)
        # stolen from https://stackoverflow.com/a/1679702
        idmap = dict(map(reversed, enumerate(ids)))
        backing = np.ndarray(n*(n-1)//2, dtype=np.float32)

        # calculate pairwise distances
        for i in range(len(ids)):
            for j in range(i+1, len(ids)):
                backing[DistMat.index(i, j, n)] = distfun(m.alns[ids[i]], m.alns[ids[j]])

        return cls(n, idmap, backing)


if __name__ == "__main__":
    import sys
    ref = [x for x in sys.argv[1]]
    alt = [x for x in sys.argv[2]]
    print(alignment_distance(ref, alt))
