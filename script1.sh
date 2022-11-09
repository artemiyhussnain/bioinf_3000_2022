#!/bin/bash

echo Running script 1...
ROOT="Assignment6"
mkdir -p ${ROOT}
cd ${ROOT}
mkdir -p data DB results
cd results
mkdir -p 1_QC 2_clean_data 3_denovo_assembly 4_genome_guided_assembly 5_final_assembly
cd ..
cd DB

echo Downloading genome sequence...
wget -O Athaliana_TAIR10.fa.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
echo Done
echo gunzipping...
gunzip Athaliana_TAIR10.fa.gz
echo Done

echo Downloading genome annotation...
wget -O Athaliana_TAIR10.gff3.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.51.gff3.gz
echo Done
echo gunzipping...
gunzip Athaliana_TAIR10.gff3.gz
echo Done

echo Downloading BUSCO lineage...
wget -O brassicales_odb10.tar.gz https://busco-data.ezlab.org/v5/data/lineages/brassicales_odb10.2020-08-05.tar.gz
echo Unzipping...
tar -zxf brassicales_odb10.tar.gz
rm brassicales_odb10.tar.gz
echo Done

echo Building genome indices for STAR and gmap...
mkdir -p Athaliana_TAIR10_STAR
STAR --runMode genomeGenerate --genomeSAindexNbases 12 --genomeDir Athaliana_TAIR10_STAR \
--genomeFastaFiles Athaliana_TAIR10.fa \
--sjdbOverhang 124 --sjdbGTFfile Athaliana_TAIR10.gff3

gmap_build -D ./ -d Athaliana_TAIR10_gmap Athaliana_TAIR10.fa
echo Done

cd ..