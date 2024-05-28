#!/bin/sh

bn=$(basename $1 .fa)

# piping seems to fail on long jobs on the cluster, so write the output to a file
~/msa-eval/carlip/msa $1 | ~/msa-eval/ultramsatric/scripts/msatofasta.py $1 > ~/msa-eval/diversealn/$bn.opt.aln.fa
