#!/bin/sh
#SBATCH --job-name cleanTest
#SBATCH -t 24:00:00
#SBATCH -p defq
#SBATCH --mail-user=ccoughlin@perimeterinstitute.ca
#SBATCH --mail-type=ALL

module load python/3.8
unset PYTHONNOUSERSITE
python testStatistics.py
