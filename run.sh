#!/bin/bash

#SBATCH --job-name="Drake_test"
#SBATCH --output="/home/sha59/4ka.out"
#SBATCH -n 64
#SBATCH --time=1000:00
#SBATCH --mem=64000
#SBATCH --mail-type=ALL

module load R
srun Rscript make.R
