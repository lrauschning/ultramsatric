#!/bin/sh

for x in blmb lyase_1 p450 rvp seatoxin; do
		./submit.sh ./regressive_call.sh $x
done
