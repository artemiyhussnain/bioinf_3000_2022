#!/bin/bash

echo Running script 2...
ROOT="Assignment6"
cd ${ROOT}
cp ~/data/Transcriptomics_data/Assignment/Col_leaf_chr4_R1.fastq.gz data
cp ~/data/Transcriptomics_data/Assignment/Col_leaf_chr4_R2.fastq.gz data

echo Getting QC for raw reads...
fastqc -t 2 -o results/1_QC data/Col_leaf_chr4_R*.fastq.gz
echo Done

echo Trimming adaptors and low quality bases...
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o results/2_clean_data/Col_leaf_chr4_R1_clean.fastq.gz -p results/2_clean_data/Col_leaf_chr4_R2_clean.fastq.gz \
-m 25 -q 20,20 data/Col_leaf_chr4_R1.fastq.gz data/Col_leaf_chr4_R2.fastq.gz
echo Done

echo Getting QC for clean reads...
cd results
fastqc -t 2 -o 1_QC 2_clean_data/Col_leaf_chr4_R*_clean.fastq.gz
echo Done