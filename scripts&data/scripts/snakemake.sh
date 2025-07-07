#!/bin/bash

#SBATCH --job-name=rna_snake
#SBATCH --mem=60G
#SBATCH --ntasks=8
#SBATCH -o /group/rbaygrp/armstrong/doseexposure/logs/out-%A.%a_snake.txt
#SBATCH -e /group/rbaygrp/armstrong/doseexposure/logs/error-%A.%a_snake.txt
#SBATCH --time=172:00:00
#SBATCH --mail-user=mlarmstrong@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -p high2

module load conda/base/latest
mamba activate snakemake

module load fastqc/0.11.9
module load salmon/1.10.1

cd /group/rbaygrp/armstrong/doseexposure/snakemake_dir

#Run this line of code if your samples were sequenced on one lane:
snakemake --latency-wait 300  --snakefile snakefile --printshellcmd --jobs 5 --conda-frontend conda

