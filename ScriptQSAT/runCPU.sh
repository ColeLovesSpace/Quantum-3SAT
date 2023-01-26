#!/bin/sh
#SBATCH --job-name plotSATscalingCPU13
#SBATCH -t 24:00:00
#SBATCH -p defq
#SBATCH -N 1
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=ccoughlin@perimeterinstitute.ca
#SBATCH --mail-type=ALL

module load python/3.8
unset PYTHONNOUSERSITE
python scalingCPU.py
