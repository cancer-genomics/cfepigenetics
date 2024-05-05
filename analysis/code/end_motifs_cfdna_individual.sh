#!/bin/bash

# Location of job execution
#$ -cwd

# Standard error stream will be merged with stdio
#$ -j y

# Job resource option: file size limit
#$ -l h_fsize=20G

# Job resource option: memory needed
#$ -l mem_free=20G

# Job resource option: max memory used
#$ -l h_vmem=20G

# Job resource option: max runtime
#$ -l h_rt=96:00:00

# Sets SGE environment: smp (shared) vs. local (in the same node), n = # cores
#$ -pe local 1

# Settings for array jobs = NUMBER OF SAMPLES
#$ -t 1-543
#$ -o /dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/cfdna_methylation_data/logs

# Max concurrent array job task execution #####
#$ -tc 200

# Choose nodes #####
#$ -l cancergen
umask g+w
module load conda_R

Rscript end_motifs_cfdna_individual.R $SGE_TASK_ID
