#!/bin/sh

mamba init
mamba activate tcoffee

export DUMP_SEQ_BUCKETS_ONLY=1
for x in blmb lyase_1 p450 rvp; do
      # output in the right dir
      mkdir $x
      cd $x

      # get the 10 most diverse sequences as the first seqdump
      t_coffee -reg -thread 4 -seq ~/homfam/combinedSeqs/${x}.fa -nseq 10

		mv $x/seqdump.1 ./${x}.div.fa

      cd ..
done

export DUMP_SEQ_BUCKETS_ONLY=0
for x in ./*.div.fa; do
      # align them using tcoffee
		t_coffee -thread 4 -seq $x -output fasta_aln -outfile $(basename $x .fa).tcof.fa
		#TODO align using msa
		# ~/msa -i ...
done
