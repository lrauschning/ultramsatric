#!/bin/sh

for x in ./*.div.fa; do
	echo Submitting ./msacall.sh $x
	./submit.sh ~/ultramsatric/scripts/msacall.sh $x
done
