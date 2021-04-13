#!/bin/bash
#
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=user@.edu
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J star_alignment
#SBATCH --output=run_test.out
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.
#SBATCH --mem=32g # Memory pool for all cores (see also --mem-per-cpu)
##array should start from zero

#adapted from https://github.com/griffithlab/rnaseq_tutorial/wiki/RNAseq-Data
#1. Read data acquisition
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
md5sum HBR_UHR_ERCC_ds_5pc.tar

#2. Read data setting
tar -xvf HBR_UHR_ERCC_ds_5pc.tar

#3. adapter trimming using cutadapt
#adapter sequence from https://cutadapt.readthedocs.io/en/stable/guide.html#paired-end
module load cutadapt/2.10
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-m 10 \  # min read length 
#-q 20 \  # min average quality per read
-o trimmed.R1.fastq.gz -p trimmed.R2.fastq.gz HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

#4. QC using fastqc
module load fastqc/0.11.9
mkdir QC
fastqc --outdir ./QC --noextract trimmed.R1.fastq.gz trimmed.R2.fastq.gz

#5. aligning reads - genomeIndex generation
module load fastqc/0.11.7
module load star/2.7.5a
module load gcc/9.2.0

STAR \
--runMode genomeGenerate \
--genomeDir ./STAR_reference/ \
--genomeFastaFiles /ihome/hpark/hyp15/test/RNASeq/ref/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile /ihome/hpark/hyp15/test/RNASeq/ref/gencode.v37.primary_assembly.annotation.gtf 

#5. aligning reads - 
STAR \
--runMode alignReads \
--outTmpDir ${FILELOC}tmp \
--runThreadN $NSLOTS \
--outSAMstrandField intronMotif \
--twopassMode Basic \
--readFilesIn trimmed.R1.fastq.gz trimmed.R2.fastq.gz \
--outFileNamePrefix $FILELOC \
--quantMode GeneCounts \
--outStd Log \
--outWigType bedGraph \
--outWigReferencesPrefix ./$FILELOC/bedgraph \
--genomeDir ./STAR_reference/ \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate
