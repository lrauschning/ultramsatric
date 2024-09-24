#!/usr/bin/env python3

import sys
import argparse


"""
Small script that takes the output of MSA (Carillo & Lipman) as stdin
and the corresponding input FASTA as argument and emits the computed
optimal MSA in FASTA format to stdout.

!!!
This script assumes the output is in one line per sample, which is not
the default behaviour of MSA.
You can make MSA output in this format by setting `LINE` in `defs.h` to
a value larger than your input sequences prior to compilation.
!!!

The input FASTA is necessary, as MSA does not preserve sample names in
the output, but does preserve sample ordering.
The diagnostic output of MSA is printed to stderr.
"""

def parse_msa(fin, fasta, fout):
    ## extract sample names from fasta
    samples = list(filter(lambda x: x[0] == '>', fasta.readlines()))

    # skip until the
    #          ***  Optimal Multiple Alignment  ***
    # line
    for line in fin:
        if "Optimal" in line:
            break

    next(fin) # skip to sequences

    # iterate through as many lines as there are samples
    for sample, seq in zip(samples, fin):
        # sample and seq both already end with newlines
        print(sample, file=fout, end='')
        print(seq, file=fout, end='')

    # now there should be an empty line
    if next(fin).strip() != '':
        print("Non-Empty line after alignment. Check that the FASTA file is correct!", file=sys.stderr)
        raise ValueError("Non-Empty line after alignment.")

    # print start of debug info
    for _ in range(4):
        print(next(fin), file=sys.stderr, end='')

    # end
    # in fin, there should only be the pairwise projected costs left over.

    



if __name__ == '__main__':
    parse_msa(sys.stdin, open(sys.argv[1], 'r'), sys.stdout)
