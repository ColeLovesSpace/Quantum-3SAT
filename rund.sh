#!/bin/sh
#SBATCH --job-name statSAT
#SBATCH -t 01:00:00
#SBATCH -p debugq
#SBATCH --mail-user=ccoughlin@perimeterinstitute.ca
#SBATCH --mail-type=ALL
module load python/3.8
unset PYTHONNOUSERSITE
python StatSATtest.py
