#!/bin/bash

echo Running script 3...
ROOT="Assignment6"
cd ${ROOT}

echo Running denovo assembly with Trinity...
mkdir -p results/3_denovo_assembly/Athaliana_trinity
Trinity --seqType fq --max_memory 8G --left results/2_clean_data/Col_leaf_chr4_R1_clean.fastq.gz \
--right results/2_clean_data/Col_leaf_chr4_R2_clean.fastq.gz --output results/3_denovo_assembly/Athaliana_trinity \
--CPU 2 --bypass_java_version_check
echo Done

echo Aligning denovo assembly to genome...
cp -r DB/Athaliana_TAIR10_gmap results/3_denovo_assembly
cd results/3_denovo_assembly
gmap -D ./ -d Athaliana_TAIR10_gmap -t 2 -f 3 -n 1 \
Athaliana_trinity/Trinity.fasta > Athaliana_trinity_gmap.gff3
rm -r Athaliana_TAIR10_gmap
echo Done

echo Assessing denovo assembly quality with BUSCO...
cd ${ROOT}
busco -f --cpu 2 -i results/3_denovo_assembly/Athaliana_trinity/Trinity.fasta -l DB/brassicales_odb10 \
--out_path results/3_denovo_assembly -o Athaliana_busco_denovo -m tran
echo Done

cd ../..