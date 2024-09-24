#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=l.rauschning@campus.lmu.de
#SBATCH --output=logs/msa-mutated-%j.out
#SBATCH --error=logs/msa-mutated-%j.err
#SBATCH --job-name=CARLIP_MSA_MUT
#SBATCH -M biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

bn=$(basename $1 .fa)

~/msa-eval/carlip/msa -c ~/msa-eval/mutated.scores $1 | ~/msa-eval/ultramsatric/scripts/msatofasta.py $1 > ~/msa-eval/diversealn/$bn.mut.aln.fa
