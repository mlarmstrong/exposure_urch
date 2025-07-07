#! /bin/bash
#SBATCH --mail-user=mlarmstrong@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=QC
#SBATCH --mem=60G
#SBATCH --partition=high2
#SBATCH --time=3-1:00:00
set -e
set -x
# To run script type:  sbatch 1.fastQC.sh

## Change directories to where the fastq files are located
cd /group/rbaygrp/armstrong/doseexposure/Data/Project_RBMA_TAG0282
# Call fastp package
module load fastqc

#RUN FOR IT MARTY
fastqc -t 2 *R1.fastq.gz -o /group/rbaygrp/armstrong/doseexposure/qual
