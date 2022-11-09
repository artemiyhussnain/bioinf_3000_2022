## Techincal note
Please use the appropriate conda environment!

`conda create -n transcriptomics_analysis`
`conda activate transcriptomics_analysis`
`conda install -c bioconda trinity busco gffread star`

## Theoretical question
*Can we get all genes assembled from one transcriptome sequenced from a sample collected from one particular tissue at one particular developmental stage?*

- No, a transcriptome assembled from sequencing data of a particular tissue will not contain all genes in the organism's genome
- This is because the pattern of expression (i.e., the few thousand genes out of the entire genome that are actually expressed) is what gives a tissue type its identity (i.e., phenotype)
- Since a transcriptome is assembled form transcribed RNA sequences, both protein-coding and non-coding, it will only contain these expressed genes