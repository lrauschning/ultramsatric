#!/bin/sh

for x in ./*.div.fa; do
	./submit.sh ./msacall.sh $x
done
