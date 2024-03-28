#!/bin/sh

cd diversealn

bn=$(basename $1 .fa)

# piping seems to fail on long jobs on the cluster, so write the output to a file
~/carlip/msa $1 > ~/diversealn/$bn.opt.out
cat ~/diversealn/$bn.opt.out | ~/ultramsatric/scripts/msatofasta.py $1 > ~/diversealn/$bn.opt.aln.fa
