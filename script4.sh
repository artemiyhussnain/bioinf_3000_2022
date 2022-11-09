#!/bin/bash

echo Running script 4...
ROOT="Assignment6"
cd ${ROOT}

echo Mapping transcripts to genome with STAR...
mkdir -p results/4_genome_guided_assembly/Athaliana_star
STAR --genomeDir DB/Athaliana_TAIR10_STAR --readFilesIn results/2_clean_data/Col_leaf_chr4_R1_clean.fastq.gz \
results/2_clean_data/Col_leaf_chr4_R2_clean.fastq.gz --readFilesCommand zcat --runThreadN 2 \
--outSAMstrandField intronMotif --outSAMattributes All --outFilterMismatchNoverLmax 0.03 --alignIntronMax 10000 \
--outSAMtype BAM SortedByCoordinate --outFileNamePrefix Col_leaf_chr4. --quantMode GeneCounts
mv Col_leaf_chr4.* results/4_genome_guided_assembly/Athaliana_star
echo Done

echo Running genome-guided assembly with StringTie...
cd results/4_genome_guided_assembly
stringtie Athaliana_star/Col_leaf_chr4.Aligned.sortedByCoord.out.bam -o Athaliana_stringtie.gtf -p 2
echo Done
cd ../..

echo Getting assembly fasta with gffread...
gffread -w results/4_genome_guided_assembly/Athaliana_stringtie.fasta -g DB/Athaliana_TAIR10.fa \
results/4_genome_guided_assembly/Athaliana_stringtie.gtf
echo Done

echo Assessing genome-guided assembly quality with BUSCO...
cd results/4_genome_guided_assembly
busco -f --cpu 2 -i Athaliana_stringtie.fasta -l ../../DB/brassicales_odb10 \
-o Athaliana_busco_genome -m tran
echo Done

echo Summarising results...
ROOT="Assignment6"
cd ${ROOT}

cp results/1_QC/Col_leaf_chr4_R1_fastqc.html results/5_final_assembly
cp results/1_QC/Col_leaf_chr4_R2_fastqc.html results/5_final_assembly
cp results/1_QC/Col_leaf_chr4_R1_clean_fastqc.html results/5_final_assembly
cp results/1_QC/Col_leaf_chr4_R2_clean_fastqc.html results/5_final_assembly
cp results/3_denovo_assembly/Athaliana_trinity/Trinity.fasta results/5_final_assembly
cp results/3_denovo_assembly/Athaliana_trinity_gmap.gff3 results/5_final_assembly
cp results/4_genome_guided_assembly/Athaliana_stringtie.fasta results/5_final_assembly
cp results/4_genome_guided_assembly/Athaliana_stringtie.gtf results/5_final_assembly

mkdir -p scripts
cd ..
cp README.md ${ROOT}
cp *.sh ${ROOT}/scripts
echo Done