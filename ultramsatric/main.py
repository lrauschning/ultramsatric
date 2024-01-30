#!/bin/env python3
#from . import __version__
import sys

from .ultrametric import *
from .msa import MSA
from .distance import *

def main():
    #print(__version__)
    print(sys.argv)
    print("===full MSA===")
    m = MSA.from_file(sys.argv[1])
    eval_msa(m)

    print("===without outgroups===")
    toxin_set = set(["1kbt","2crt","2cdx","1cdta","1tgxa","1kxia","1tfs","1drs","1txb","2ctx","1ntn","1lsi","2abxa","2nbta","1nean","1nor","1cod","1nxb","1ntx","1fas"])
    m = m.subset(toxin_set)
    eval_msa(m)

    #print("===seatoxin subset===")
    #toxin_set = set(["1kbt","2crt","2cdx","1cdta","1tgxa","1kxia","1tfs","1drs","1txb","2ctx","1ntn","1lsi","2abxa","2nbta","1nean","1nor","1cod","1nxb","1ntx","1fas"])
    #m = m.subset(toxin_set)
    #eval_msa(m)

def eval_msa(m: MSA):
    #print(m)
    d = DistMat.from_msa(m, distfun=scoredist)
    #d = DistMat.from_msa(m, distfun=alignment_distance)
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
    print("SOP", "", "sum:", "",  m.sumofpairs(scoredist), sep='\t')
    print("SOP", "", "avg:", "", m.sumofpairs_avg(scoredist), sep='\t')



if __name__ == "__main__":
    main()
