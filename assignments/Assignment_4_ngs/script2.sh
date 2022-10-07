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

if [[ -f ${RAW_FASTQ}/${NAME}_1.fastq.gz ]]; then
echo Sequencing data already downloaded
else
echo -e Downloading sequencing data to ${RAW_FASTQ}'...'
wget -O ${RAW_FASTQ}/${NAME}_1.fastq.gz \
https://universityofadelaide.box.com/shared/static/egl3n16r0ziaxlvbs9074xqd1liktnuz.gz
wget -O ${RAW_FASTQ}/${NAME}_2.fastq.gz \
https://universityofadelaide.box.com/shared/static/g2ly4kzz1blus5juy426i37zl45o38pu.gz
echo Done
fi

echo -e '\nGenerating raw read FASTQC report at '${RAW_FASTQC}'...'
fastqc -t 2 -o ${RAW_FASTQC} ${RAW_FASTQ}/*.gz
echo Done

Bases with PHRED score below 28 are trimmed from both 5' and 3' ends
Illumina TruSeq paired-read adapter sequences can be found at \
https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
echo -e '\nSaving reads with trimmed poor quality bases and adapters to '${TRIMMED_FASTQ}'...'
cutadapt \
-m 35 \
-q 28,28 \
-a CTGTCTCTTATACACATCT \
-A CTGTCTCTTATACACATCT \
-o ${TRIMMED_FASTQ}/${NAME}_1-trimmed.fastq.gz \
-p ${TRIMMED_FASTQ}/${NAME}_2-trimmed.fastq.gz \
${RAW_FASTQ}/${NAME}_1.fastq.gz ${RAW_FASTQ}/${NAME}_2.fastq.gz > \
${TRIMMED_LOG}/cutadapt.log
echo Done

echo -e '\nGenerating trimmed read FASTQC report at '${TRIMMED_FASTQ}'...'
fastqc -t 2 -o ${TRIMMED_FASTQC} ${TRIMMED_FASTQ}/*.gz
echo Done

#Excluding unmapped reads (SAM flag 4)
echo -e '\nSaving trimmed read and genome alignment to '${ALIGNED}'...'
bwa mem -t 2 ${REFS_ATHALIANA}/bwa_index/Athaliana.TAIR10-51_index ${TRIMMED_FASTQ}/${NAME}_1-trimmed.fastq.gz \
${TRIMMED_FASTQ}/${NAME}_2-trimmed.fastq.gz | \
samtools view -bhS -F4 - > ${ALIGNED}/${NAME}-aligned.bam
echo Done

echo -e '\nSorting and indexing alignment at '${ALIGNED}'...'
samtools sort -o ${ALIGNED}/${NAME}-aligned.bam ${ALIGNED}/${NAME}-aligned.bam
samtools index ${ALIGNED}/${NAME}-aligned.bam
echo Done