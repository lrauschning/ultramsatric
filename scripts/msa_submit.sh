#!/bin/sh

for x in ./*.div.fa; do
	echo Submitting ./msa_call.sh $x
	#./submit.sh ~/ultramsatric/scripts/msa_call.sh $x
	echo Submitting ./msa_mutatedcall.sh $x
	./submit.sh ~/ultramsatric/scripts/msa_mutatedcall.sh $x
done
