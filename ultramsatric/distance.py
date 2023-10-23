#!/bin/python3
from typing import List, Callable
import numpy as np


def identity(ref:chr, alt:chr) -> float:
    return 1 if ref != alt else 0

def blosum(ref:chr, alt:chr) -> float:
    pass

def linear(n:int) -> float:
    return 3*n

def affine(n:int) -> float:
    return 3 + 2*n

def alignment_distance(ref: List[chr], alt: List[chr], dfun:Callable[[chr, chr], float] = identity, gapcost:Callable[[int], float] = linear) -> float:
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


#blosum = #blosum matrix as np.ndarray

if __name__ == "__main__":
    import sys
    ref = [x for x in sys.argv[1]]
    alt = [x for x in sys.argv[2]]
    print(alignment_distance(ref, alt))
