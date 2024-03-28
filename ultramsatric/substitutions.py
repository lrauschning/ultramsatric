#!/bin/python3

from typing import List, Callable, Dict
import math
import itertools

from collections import defaultdict

import blosum as bl

BLOSUM = bl.BLOSUM(62)

PAM250 = {
        'A': {'A':2, 'R':-2, 'N':0, 'D':0, 'C':-2, 'Q':0, 'E':0, 'G':1, 'H':-1, 'I':-1, 'L':-2, 'K':-1, 'M':-1, 'F':-3, 'P':1, 'S':1, 'T':1, 'W':-6, 'Y':-3, 'V':0},
        'R': {'A':-2, 'R':6, 'N':0, 'D':-1, 'C':-4, 'Q':1, 'E':-1, 'G':-3, 'H':2, 'I':-2, 'L':-3, 'K':3, 'M':0, 'F':-4, 'P':0, 'S':0, 'T':-1, 'W':2, 'Y':-4, 'V':-2},
        'N': {'A':0, 'R':0, 'N':2, 'D':2, 'C':-4, 'Q':1, 'E':1, 'G':0, 'H':2, 'I':-2, 'L':-3, 'K':1, 'M':-2, 'F':-3, 'P':0, 'S':1, 'T':0, 'W':-4, 'Y':-2, 'V':-2},
        'D': {'A':0, 'R':-1, 'N':2, 'D':4, 'C':-5, 'Q':2, 'E':3, 'G':1, 'H':1, 'I':-2, 'L':-4, 'K':0, 'M':-3, 'F':-6, 'P':-1, 'S':0, 'T':0, 'W':-7, 'Y':-4, 'V':-2},
        'C': {'A':-2, 'R':-4, 'N':-4, 'D':-5, 'C':12, 'Q':-5, 'E':-5, 'G':-3, 'H':-3, 'I':-2, 'L':-6, 'K':-5, 'M':-5, 'F':-4, 'P':-3, 'S':0, 'T':-2, 'W':-8, 'Y':0, 'V':-2},
        'Q': {'A':0, 'R':1, 'N':1, 'D':2, 'C':-5, 'Q':4, 'E':2, 'G':-1, 'H':3, 'I':-2, 'L':-2, 'K':1, 'M':-1, 'F':-5, 'P':0, 'S':-1, 'T':-1, 'W':-5, 'Y':-4, 'V':-2},
        'E': {'A':0, 'R':-1, 'N':1, 'D':3, 'C':-5, 'Q':2, 'E':4, 'G':0, 'H':1, 'I':-2, 'L':-3, 'K':0, 'M':-2, 'F':-5, 'P':-1, 'S':0, 'T':0, 'W':-7, 'Y':-4, 'V':-2},
        'G': {'A':1, 'R':-3, 'N':0, 'D':1, 'C':-3, 'Q':-1, 'E':0, 'G':5, 'H':-2, 'I':-3, 'L':-4, 'K':-2, 'M':-3, 'F':-5, 'P':0, 'S':1, 'T':0, 'W':-7, 'Y':-5, 'V':-1},
        'H': {'A':-1, 'R':2, 'N':2, 'D':1, 'C':-3, 'Q':3, 'E':1, 'G':-2, 'H':6, 'I':-2, 'L':-2, 'K':0, 'M':-2, 'F':-2, 'P':0, 'S':-1, 'T':-1, 'W':-3, 'Y':0, 'V':-2},
        'I': {'A':-1, 'R':-2, 'N':-2, 'D':-2, 'C':-2, 'Q':-2, 'E':-2, 'G':-3, 'H':-2, 'I':5, 'L':2, 'K':-2, 'M':2, 'F':1, 'P':-2, 'S':-1, 'T':0, 'W':-5, 'Y':-1, 'V':4},
        'L': {'A':-2, 'R':-3, 'N':-3, 'D':-4, 'C':-6, 'Q':-2, 'E':-3, 'G':-4, 'H':-2, 'I':2, 'L':6, 'K':-3, 'M':4, 'F':2, 'P':-3, 'S':-3, 'T':-2, 'W':-2, 'Y':-1, 'V':2},
        'K': {'A':-1, 'R':3, 'N':1, 'D':0, 'C':-5, 'Q':1, 'E':0, 'G':-2, 'H':0, 'I':-2, 'L':-3, 'K':5, 'M':0, 'F':-5, 'P':-1, 'S':0, 'T':0, 'W':-3, 'Y':-4, 'V':-2},
        'M': {'A':-1, 'R':0, 'N':-2, 'D':-3, 'C':-5, 'Q':-1, 'E':-2, 'G':-3, 'H':-2, 'I':2, 'L':4, 'K':0, 'M':6, 'F':0, 'P':-2, 'S':-2, 'T':-1, 'W':-4, 'Y':-2, 'V':2},
        'F': {'A':-3, 'R':-4, 'N':-3, 'D':-6, 'C':-4, 'Q':-5, 'E':-5, 'G':-5, 'H':-2, 'I':1, 'L':2, 'K':-5, 'M':0, 'F':9, 'P':-5, 'S':-3, 'T':-3, 'W':0, 'Y':7, 'V':-1},
        'P': {'A':1, 'R':0, 'N':0, 'D':-1, 'C':-3, 'Q':0, 'E':-1, 'G':0, 'H':0, 'I':-2, 'L':-3, 'K':-1, 'M':-2, 'F':-5, 'P':6, 'S':1, 'T':0, 'W':-6, 'Y':-5, 'V':-1},
        'S': {'A':1, 'R':0, 'N':1, 'D':0, 'C':0, 'Q':-1, 'E':0, 'G':1, 'H':-1, 'I':-1, 'L':-3, 'K':0, 'M':-2, 'F':-3, 'P':1, 'S':2, 'T':1, 'W':-2, 'Y':-3, 'V':-1},
        'T': {'A':1, 'R':-1, 'N':0, 'D':0, 'C':-2, 'Q':-1, 'E':0, 'G':0, 'H':-1, 'I':0, 'L':-2, 'K':0, 'M':-1, 'F':-3, 'P':0, 'S':1, 'T':3, 'W':-5, 'Y':-3, 'V':0},
        'W': {'A':-6, 'R':2, 'N':-4, 'D':-7, 'C':-8, 'Q':-5, 'E':-7, 'G':-7, 'H':-3, 'I':-5, 'L':-2, 'K':-3, 'M':-4, 'F':0, 'P':-6, 'S':-2, 'T':-5, 'W':17, 'Y':0, 'V':-6},
        'Y': {'A':-3, 'R':-4, 'N':-2, 'D':-4, 'C':0, 'Q':-4, 'E':-4, 'G':-5, 'H':0, 'I':-1, 'L':-1, 'K':-4, 'M':-2, 'F':7, 'P':-5, 'S':-3, 'T':-3, 'W':0, 'Y':10, 'V':-2},
        'V': {'A':0, 'R':-2, 'N':-2, 'D':-2, 'C':-2, 'Q':-2, 'E':-2, 'G':-1, 'H':-2, 'I':4, 'L':2, 'K':-2, 'M':2, 'F':-1, 'P':-1, 'S':-1, 'T':0, 'W':-6, 'Y':-2, 'V':4},
        }
# from https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/PAM250

# generated using the following python code:
# f = open('pam.tsv') 
# x = f.readline()
# ind = x.strip().split('\t')
# for i, l in enumerate(f):
#   print("'", ind[i], "'", ': {', ', '.join(["'" + ind[i] + "'" + ':' + v for i, v in enumerate(l.strip().split('\t')[1:])]), '},', sep='')

AA_FREQS = {'A': 8.76, 'R': 5.78, 'N': 3.93, 'D': 5.49, 'C': 1.38, 'Q': 3.9,
           'E': 6.32, 'G': 7.03, 'H': 2.26, 'I': 5.49, 'L': 9.68, 'K': 5.19, 'M': 2.32,
           'F': 3.87, 'P': 5.02, 'S': 7.14, 'T': 5.53, 'W': 1.25, 'Y': 2.91, 'V': 6.73
        }
# from https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties

## Substitution scores

def identity(ref:chr, alt:chr) -> float:
    return 1 if ref != alt else 0

def blosum(ref:chr, alt:chr) -> float:
    return BLOSUM[ref][alt]

def pam(ref:chr, alt:chr) -> float:
    return PAM250[ref][alt]

def from_msa_format(fin):
    """
    Takes a substitution model in the format used by MSA (Carillo & Lipman),
    returns a function corresponding to it. Assumes the input to be symmetric, like MSA.
    Ignores the gapcost specified in the format.
    Empty lines and lines starting with '#' after the first line are ignored.
    """
    lookup = defaultdict(lambda: dict())
    next(fin)
    for line in fin:
        if line.strip() == '' or line[0] == '#': # ignore empty and comment lines
            continue

        fields = line.strip().split(' ')
        assert(len(fields) == 3) # check that the format is correct

        if fields[0] <= fields[1]: # enforce lexicalic sorting in the dictionary
            lookup[fields[0]][fields[1]] = float(fields[2])
        else: 
            lookup[fields[1]][fields[0]] = float(fields[2])

    return lambda x, y: lookup[x][y] if x <= y else lookup[y][x]


def get_ev(subs: Callable[[chr,chr], float], eqdist: bool = False) -> float:
    """
    Computes the per-position expectation value of a substitution model.
    Uses AA frequencies from a database by default, if eqdist is set to `True` assumes an even distribution of AAs.
    """
    #aas = [x for x in BLOSUM][:-5] # remove values coding for unknown AAs
    freqs = itertools.product(AA_FREQS.items(), AA_FREQS.items()) if not eqdist\
            else [1/len(AA_FREQS)**2 for _ in AA_FREQS]
    return sum(map(lambda x: x[0][1]*x[1][1]*subs(x[0][0], x[1][0]),
        freqs)) / (len(AA_FREQS)**2)

## Gapcost functions
def linear_cons(m:float) -> Callable[int, float]:
    return lambda n: m*n

def linear(n:int) -> float:
    return 3*n

def affine_cons(m:float, t:float) -> Callable[int, float]:
    return lambda n: m*n + t

def affine(n:int) -> float:
    return 3 + 2*n

def no_gaps_cons() -> Callable[int, float]:
    return lambda n: 0.0

def no_gaps(n:int) -> float:
    return 0.0

