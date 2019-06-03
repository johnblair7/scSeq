#!/bin/bash
## SLURM parameters
#SBATCH --time 3-0  # days-minutes
#SBATCH --job-name=d120_all_seurat # Job name
#SBATCH --account=fc_jbscseq
#SBATCH --nodes=1
#SBATCH --ntasks=19 # Number of cores
#SBATCH --mem=64000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=savio # Partition to submit to
#SBATCH --qos=savio_normal
#SBATCH --output=jobs/countArrayJob_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=jobs/countArrayJob_%A_%a.err # File to which STDERR will be written

module load r/3.5.1


NOW=$(date +"_%m%d%Y-%H%M%S")

R_LIBS=/global/scratch/johnblair/rlibs/3.5 R --vanilla < SeuratClusterManual24.R  > 'Seurat120'$NOW'.Rout'
