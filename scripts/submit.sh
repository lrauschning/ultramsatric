#!/bin/sh

qsub -pe smp 10 -l virtual_free=64G -q cn-el7 tcoffee.sh $1