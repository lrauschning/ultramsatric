#!/bin/sh

cd diversealn

~/carlip/msa -c ~/ultramsatric/scripts/mutated.scores $1 | ~/ultramsatric/scripts/msatofasta.py $1 > ~/diversealn/$(basename $1 .fa).mut.aln
