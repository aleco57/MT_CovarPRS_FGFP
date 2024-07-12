#!/bin/bash
#SBATCH --job-name="scripts01_02"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1 
#SBATCH --time=02:00:00
#SBATCH --mem=20G
#SBATCH --account=sscm009461

module add languages/r/4.3.1

Rscript 01_ExtractPhenoIns.R

source 02_PRSice.sh  
