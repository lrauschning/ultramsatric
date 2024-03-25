#!/bin/sh

mamba activate tcoffee

export DUMP_SEQ_BUCKETS_ONLY=1

#cd diversealn

# first two added manually,
# rest obtained using
# ls ../homfam/combinedSeqs | shuf | head -n 15
cat ../homfam/combinedSeqs/seatoxin.fa \
		../homfam/combinedSeqs/adh.fa \
		../homfam/combinedSeqs/subt.fa \
		../homfam/combinedSeqs/sodcu.fa \
		../homfam/combinedSeqs/toxin.fa \
		../homfam/combinedSeqs/myb_DNA-binding.fa \
		../homfam/combinedSeqs/Acetyltransf.fa \
		../homfam/combinedSeqs/annexin.fa \
		../homfam/combinedSeqs/blmb.fa \
		../homfam/combinedSeqs/seatoxin.fa \
		../homfam/combinedSeqs/DMRL_synthase.fa \
		../homfam/combinedSeqs/ghf10.fa \
		../homfam/combinedSeqs/az.fa \
		../homfam/combinedSeqs/ace.fa \
		../homfam/combinedSeqs/hip.fa \
		../homfam/combinedSeqs/lyase_1.fa \
| t_coffee -reg -thread 10 -seq STDIN -nseq 50 -outfile alignment.aln
