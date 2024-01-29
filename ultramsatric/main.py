#!/bin/env python3
#from . import __version__
import sys

from ultramsatric import *
from .ultrametric import *
from .msa import MSA
from .distance import *

def main():
    #print(__version__)

    print(sys.argv)
    m = MSA.from_file(sys.argv[1])
    #print(m)
    d = DistMat.from_msa(m, distfun=scoredist)
    #print(d.idmap)
    #print(d)
    d_upgma = UPGMA_matrix(d)
    #print(d_upgma.idmap)
    #print(d_upgma)
    diff = d - d_upgma
    #print(diff)
    print("msatric", "frobenius:", diff.norm_frobenius(), sep='\t')
    print("msatric", "absavg:", diff.absavg(), sep='\t')
    print("distance", "frobenius:", d.norm_frobenius(), sep='\t')
    print("distance", "absavg:", d.absavg(), sep='\t')
    print("SOP", "", "sum:", "",  m.sumofpairs(), sep='\t')
    print("SOP", "", "avg:", "", m.sumofpairs_avg(), sep='\t')



if __name__ == "__main__":
    main()
