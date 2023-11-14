#!/bin/env python3
#from . import __version__
import sys

from msa import MSA
from distance import *
from ultrametric import *

def main():
    #print(__version__)

    print(sys.argv)
    m = MSA.from_file(sys.argv[1])
    #print(m)
    d = DistMat.from_msa(m, distfun=scoredist)
    print(d.idmap)
    print(d)
    d_upgma = UPGMA_matrix(d)
    print(d_upgma.idmap)
    print(d_upgma)
    diff = d - d_upgma
    print(diff)
    print(diff.norm_frobenius())
    print(diff.absavg())



if __name__ == "__main__":
    main()
