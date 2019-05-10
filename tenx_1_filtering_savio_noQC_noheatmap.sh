#!/bin/bash
## SLURM parameters
#SBATCH --time 3-0  # days-minutes
#SBATCH --job-name=d120_etoh # Job name
#SBATCH --account=fc_jbscseq
#SBATCH --nodes=1
#SBATCH --ntasks=19 # Number of cores
#SBATCH --mem=64000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=savio # Partition to submit to
#SBATCH --qos=savio_normal
#SBATCH --output=jobs/countArrayJob_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=jobs/countArrayJob_%A_%a.err # File to which STDERR will be written

export PATH=/global/home/users/johnblair/gsl/bin:$PATH
export LD_LIBRARY_PATH=/global/home/users/johnblair/gsl/lib/:$LD_LIBRARY_PATH

module load r/3.5.1

ncores=19
NOW=$(date +"_%m%d%Y-%H%M%S")

R_LIBS=/global/scratch/johnblair/rlibs/3.5 R --vanilla < tenx_1_filtering_posconheatmap.R --args --expt d120_etoh --ncores $ncores --aggr AGGR_120_etOH --exptinfo d120_etoh_expt_info.csv --annotation newref_Ai9 --hkfile /global/scratch/johnblair/scseq_analysis/ref/ai9_negcon.txt --posctrlfile /global/scratch/johnblair/scseq_analysis/ref/ai9_PosCon.txt --runQC FALSE  > 'tenx_1_filtering'$NOW'.Rout'
