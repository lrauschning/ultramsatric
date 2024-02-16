#!/bin/env python3
from . import __version__
from .ultrametric import *
from .msa import MSA
from .distance import *

import argparse as ap
import numpy as np

def main():

    parser = ap.ArgumentParser(description="""
    ultramsatric – evaluate MSAs based on their ultrametricity
    """)
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument("-i", dest='infile', default='-', type=ap.FileType('r'), help="Input MSA in FASTA format. Default stdin. TODO support gzip input.")
    parser.add_argument("-o", dest='outfile', default='-', type=ap.FileType('wt'), help="File to write output CSV to. Default stdout.")
    parser.add_argument("-m", "--metrics", dest='metrics', default='ufrob,uabsavg', type=str, help="Metrics to compute, separated by ','. The order of metrics will be preserved in the output CSV. Valid metrics are 'frob', 'absavg', 'dfrob', 'dabsavg'. Metrics starting with 'd' are run on the distance matrix directly instead of the matrix containing the distance to the closest ultrametric tree. Default 'frob,absavg'.")
    parser.add_argument("-d", "--dist", "--distance", dest='dist', default='scoredist', type=str, help="Distance function to use to calculate a distance matrix from an MSA. Default scoredist. Can be 'scoredist', 'alndist' or 'logalndist'.")
    parser.add_argument("--id", dest='id', default=None, type=str, help="Sample ID to index the CSV with")
    parser.add_argument("--no-header", dest='header', action='store_false', default=True, help="Emit a CSV without a header")
    parser.add_argument("-p", "--print-matrix", dest='print_matrix', action='store_true', default=False, help="Print the raw matrices caculated by ultramsatric.")

    args = parser.parse_args()

        
    distmapper = {'scoredist': scoredist,
                  'alndist': alndist,
                  'logalndist': log_alndist}
    if args.dist not in distmapper:
        raise ValueError(f"Invalid argument passed to -d: {args.dist}")
    else:
        args.dist = distmapper[args.dist]

    m = MSA.from_inputstream(args.infile)

    d = DistMat.from_msa(m, distfun=args.dist)

    d_upgma = UPGMA_matrix(d)
    d_nj = NJ_matrix(d)
    d_root = root_ext_add(d)
    d_tallest = tallest_ultrametric(d)


    udiff = d - d_upgma
    ndiff = d - d_nj
    rdiff = d - d_root
    tdiff = d - d_tallest

    metricmapper = {'ufrob': lambda: str(udiff.norm_frobenius()),
                    'uabsavg': lambda: str(udiff.absavg()),
                    'ucorr': lambda: str(d.corr(d_upgma)),
                    'tfrob': lambda: str(tdiff.norm_frobenius()),
                    'tabsavg': lambda: str(tdiff.absavg()),
                    'tcorr': lambda: str(d.corr(d_tallest)),
                    'rfrob': lambda: str(rdiff.norm_frobenius()),
                    'rabsavg': lambda: str(rdiff.absavg()),
                    'rcorr': lambda: str(d.corr(d_root)),
                    'nfrob': lambda: str(ndiff.norm_frobenius()),
                    'nabsavg': lambda: str(ndiff.absavg()),
                    'ncorr': lambda: str(d.corr(d_nj)),
                    'dfrob': lambda: str(d.norm_frobenius()),
                    'dabsavg': lambda: str(d.absavg())
                    }

    if args.print_matrix:
        print("===UPGMA Matrix===")
        print(udiff)
        print("===NJ Matrix===")
        print(ndiff)
        print("===Rooting Matrix===")
        print(rdiff)
        print("===Tallest Ultrametirc Matrix===")
        print(tdiff)

    if args.header:
        if args.id:
            args.outfile.write('id,')
        args.outfile.write(args.metrics)
        args.outfile.write('\n')
    if args.id:
        args.outfile.write(args.id)
        args.outfile.write(',')
    args.outfile.write(','.join([metricmapper[x.strip()]() for x in args.metrics.split(',')]))
    args.outfile.write('\n')



    #print("===without outgroups===")
    #toxin_set = set(["1kbt","2crt","2cdx","1cdta","1tgxa","1kxia","1tfs","1drs","1txb","2ctx","1ntn","1lsi","2abxa","2nbta","1nean","1nor","1cod","1nxb","1ntx","1fas"])


#def (m: MSA) -> float:
#    return get_diffmatrix(m).norm_frobenius()
    #print("msatric", "frobenius:", diff.norm_frobenius(), sep='\t')
    #print("msatric", "absavg:", diff.absavg(), sep='\t')
    #print("distance", "frobenius:", d.norm_frobenius(), sep='\t')
    #print("distance", "absavg:", d.absavg(), sep='\t')
    #print("SOP", "", "sum:", "",  m.sumofpairs(scoredist), sep='\t')
    #print("SOP", "", "avg:", "", m.sumofpairs_avg(scoredist), sep='\t')



if __name__ == "__main__":
    main()
