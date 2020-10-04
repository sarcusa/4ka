#!/bin/bash

#SBATCH --time=1000:00
#SBATCH --mem=64000
#SBATCH --mail-type=ALL

module load R
srun Rscript test_BS.R

