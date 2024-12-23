#!/bin/sh
#SBATCH --job-name=CARLIP_MSA
#SBATCH --time=1-0:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

for x in ./*.div.fa; do
	echo Submitting ./msa_call.sh $x
	sbatch ~/msa-eval/ultramsatric/scripts/msa_call.sh $x
	#./submit.sh ~/ultramsatric/scripts/msa_call.sh $x
	echo Submitting ./msa_mutatedcall.sh $x
	sbatch ~/msa-eval/ultramsatric/scripts/msa_mutatedcall.sh $x
	#./submit.sh ~/ultramsatric/scripts/msa_mutatedcall.sh $x
done
