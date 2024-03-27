#!/bin/sh

~/carlip/msa $1 | ~/ultramsatric/scripts/msatofasta.py $1 > ~/diversealn/$(basename $1 .fa).opt.aln
