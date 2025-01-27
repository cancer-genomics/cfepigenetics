#!/bin/bash
#SBATCH --job-name=slurm_out
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=96:00:00
#SBATCH --partition=cancergen,shared
#SBATCH --output=/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/cfdna_methylation_data/logs/%x.o%A.%a.txt
#SBATCH --array=1-543%200
umask g+w
module load conda_R

Rscript Loyfer_cases.R $SLURM_ARRAY_TASK_ID
