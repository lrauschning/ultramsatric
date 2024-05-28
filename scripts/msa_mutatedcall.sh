#!/bin/sh
#SBATCH --job-name=CARLIP_MSA
#SBATCH --time=1d-0:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

bn=$(basename $1 .fa)

~/msa-eval/carlip/msa -c ~/msa-eval/mutated.scores $1 | ~/msa-eval/ultramsatric/scripts/msatofasta.py $1 > ~/msa-eval/diversealn/$bn.mut.aln.fa
