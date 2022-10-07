#!/bin/bash

NAME=SRR5882797_10M
GENOME=Athaliana.TAIR10-51
ROOT=~/Assignment4
mkdir -p ${ROOT}
cd ${ROOT}
REFS_ATHALIANA=${ROOT}/Refs/Athaliana
mkdir -p ${REFS_ATHALIANA}/bwa_index
RAW=${ROOT}/01_rawData
RAW_FASTQ=${RAW}/fastq
mkdir -p ${RAW_FASTQ}
RAW_FASTQC=${RAW}/fastqc
mkdir -p ${RAW_FASTQC}
TRIMMED=${ROOT}/02_trimmedData
TRIMMED_FASTQ=${TRIMMED}/fastq
mkdir -p ${TRIMMED_FASTQ}
TRIMMED_FASTQC=${TRIMMED}/fastqc
mkdir -p ${TRIMMED_FASTQC}
TRIMMED_LOG=${TRIMMED}/log
mkdir -p ${TRIMMED_LOG}
ALIGNED=${ROOT}/03_alignedData/bam
mkdir -p ${ALIGNED}
echo Defined key directories...

if [[ -f ${REFS_ATHALIANA}/${GENOME}.fa.gz ]]; then
echo -e '\nGenome sequence already downloaded'
else
echo -e '\nDownloading genome sequence...'
wget -O ${REFS_ATHALIANA}/${GENOME}.fa.gz \
ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
echo Done
fi

if [[ -f ${REFS_ATHALIANA}/${GENOME}.gff3.gz ]]; then
echo -e '\nGenome annotation already downloaded'
else
echo -e '\nDownloading genome annotation...'
wget -O ${REFS_ATHALIANA}/${GENOME}.gff3.gz \
ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.51.gff3.gz
echo Done
fi

echo -e '\nCalculating number of chromosomes in genome...'
cp ${REFS_ATHALIANA}/${GENOME}.fa.gz ${REFS_ATHALIANA}/${GENOME}.tmp.fa.gz
gunzip ${REFS_ATHALIANA}/${GENOME}.tmp.fa.gz
#Chromosome numbers are prefaced with > in a .fa file
#The .tmp.fa file also contains excluded sequences of Mt (mitochondrial) and Pt (plastid) chromosomes
CHROMOSOME_NUM=$(cat ${REFS_ATHALIANA}/${GENOME}.tmp.fa | grep -e '^>[0-9]\s' | wc -l)
echo ${CHROMOSOME_NUM}

echo -e '\nCalculating number of unique protein-coding genes...'
cp ${REFS_ATHALIANA}/${GENOME}.gff3.gz ${REFS_ATHALIANA}/${GENOME}.tmp.gff3.gz
gunzip ${REFS_ATHALIANA}/${GENOME}.tmp.gff3.gz
#Types of annotated elements are in 3rd tab-separated field of a .gff3 file after a 16-line header
GENE_NUM=$(cat ${REFS_ATHALIANA}/${GENOME}.tmp.gff3 | cut -f3 | tail +16 | grep -e '^gene' | wc -l)
rm ${REFS_ATHALIANA}/${GENOME}.tmp.gff3
echo ${GENE_NUM}

echo -e '\nCreating genome sequence index...'
bwa index ${REFS_ATHALIANA}/${GENOME}.tmp.fa -p ${REFS_ATHALIANA}/bwa_index/${GENOME}_index
rm ${REFS_ATHALIANA}/${GENOME}.tmp.fa
echo Done