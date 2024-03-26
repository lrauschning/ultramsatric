#!/bin/sh

for x in blmb lyase_1 p450 rvp; do
      # output in the right dir
      mkdir $x
      cd $x

      # get the 50 most diverse sequences
      export DUMP_SEQ_BUCKETS_ONLY=1
      t_coffee -reg -thread 10 -seq ~/homfam/combinedSeqs/${x}.fa -nseq 50

      # align them
      export DUMP_SEQ_BUCKETS_ONLY=0
      t_coffee -thread 10 -seq ./seqdump.1 -output fasta_aln -outfile diverse.aln.fa

      cd ..
done
