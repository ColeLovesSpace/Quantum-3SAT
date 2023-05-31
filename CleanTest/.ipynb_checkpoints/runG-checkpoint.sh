#!/bin/sh
#SBATCH --job-name plotSATscalingGPU
#SBATCH -t 01:00:00
#SBATCH -p gpudebugq
#SBATCH --gpus 1
#SBATCH --mail-user=ccoughlin@perimeterinstitute.ca
#SBATCH --mail-type=ALL
module load cuda
module load python/3.8
unset PYTHONNOUSERSITE
python GPUqtest.py
