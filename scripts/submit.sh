#!/bin/sh

# ensure the job is executed in the right dir
cd ~/diversealn

qsub -pe smp 10 -l virtual_free=64G -q cn-el7 $1 $2
