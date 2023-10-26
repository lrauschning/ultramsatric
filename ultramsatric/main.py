#!/bin/env python3
#from . import __version__
import sys

from msa import MSA
from distance import *

def main():
    print("Hello world!")
    #print(__version__)

    print(sys.argv)
    m = MSA.from_file(sys.argv[1])
    print(m)
    d = DistMat.from_msa(m, distfun=alignment_distance)
    print(d)


if __name__ == "__main__":
    main()
