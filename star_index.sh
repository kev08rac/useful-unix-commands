#!/bin/bash
#
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=user@.edu
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J star_index
#SBATCH --output=run_test.out
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.
#SBATCH --mem=8g # Memory pool for all cores (see also --mem-per-cpu)

module load fastqc/0.11.7
module load star/2.7.5a
module load gcc/9.2.0

STAR \
--runMode genomeGenerate \
--genomeDir ./STAR_reference/ \
--genomeFastaFiles /ihome/hpark/hyp15/test/RNASeq/ref/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile /ihome/hpark/hyp15/test/RNASeq/ref/gencode.v37.primary_assembly.annotation.gtf 
