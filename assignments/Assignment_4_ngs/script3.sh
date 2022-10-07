#!/bin/bash

NAME=SRR5882797_10M
ROOT=~/Assignment4
REFS_ATHALIANA=${ROOT}/Refs/Athaliana
RAW=${ROOT}/01_rawData
RAW_FASTQ=${RAW}/fastq
RAW_FASTQC=${RAW}/fastqc
TRIMMED=${ROOT}/02_trimmedData
TRIMMED_FASTQ=${TRIMMED}/fastq
TRIMMED_FASTQC=${TRIMMED}/fastqc
TRIMMED_LOG=${TRIMMED}/log
ALIGNED=${ROOT}/03_alignedData/bam

echo -e '#Summary of read alignment '${ALIGNED}'/'${NAME}'-aligned.bam' > ${ROOT}/alignmentSummary.tsv
echo -e 'Saving alignment summary to '${ROOT}'...'
#Each read is a new line
READS=$(samtools view ${ALIGNED}/${NAME}-aligned.bam | wc -l)
echo -e 'total_reads\t'${READS} >> ${ROOT}/alignmentSummary.tsv
#Only including paired reads (SAM flag 1)
PAIRS=$(samtools view -f 1 ${ALIGNED}/${NAME}-aligned.bam | wc -l)
echo -e 'total_pairs\t'${PAIRS} >> ${ROOT}/alignmentSummary.tsv
#Only including properly paired reads (SAM flag 2)
PROPER_PAIRS=$(samtools view -f 2 ${ALIGNED}/${NAME}-aligned.bam | wc -l)
echo -e 'total_proper_pairs\t'${PROPER_PAIRS} >> ${ROOT}/alignmentSummary.tsv
echo Done