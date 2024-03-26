#!/bin/sh

mamba init
mamba activate tcoffee

export DUMP_SEQ_BUCKETS_ONLY=1
mkdir $1
cd $1

# get the 10 most diverse sequences as the first seqdump
t_coffee -reg -thread 10 -seq ~/homfam/combinedSeqs/${x}.fa -nseq 10

mv ./seqdump.1 ../${x}.div.fa
