#!/bin/bash
#
#SBATCH --job-name=fastqc
#SBATCH -n 1
#SBATCH -t 0-00:60 # Runtime in D-HH:MM
#SBATCH --output=fastqc.out

module load fastqc/0.11.7

fastqc -o ./fastqc_pretrim/ early/*_1.fastq.gz
fastqc -o ./fastqc_pretrim/ late/*_2.fastq.gz