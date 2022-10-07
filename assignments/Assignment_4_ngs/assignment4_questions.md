### Question 1
*What is the difference between an index and a barcode?*

##### Indices

- An index is a sequence added to the sample (short fragment of the amplified organism DNA), along with sequencing primers and sequences complementary to the flow cell oligos.
- If multiple organisms/cell lines are being sequenced, samples coming from each one will contain unique identifying indices for both forward and reverse read direction
- In the sequencing machine, index primers are added and the index regions sequenced to identify each read before sequencing of the actual sample can begin
- Finally, the sequencing primers are added and reads are generated

##### Barcodes

- While indices are identifier sequences in the adapter regions, barcodes are added in-line, directly between the sequencing primer and sample sequence proper
- Different barcodes can be added to identify both forward and reverse reads
- They will be included in the .fastq file containing reads, and will have to be identified and removed in a pipeline


### Question 2
*What is the difference between a SAM file and a BAM file?*

- A SAM file is a human-readable plain text file that stores sequencing reads aligned to a reference genome
- In addition to the read sequences, it contains tab-delimited information on read quality, pairing with reverse read, position along reference genome, and presence of indels and substitutions
- A BAM file contains the exact same information in binary
- This reduces overall file size, but is not human-readable
- SAM/BAM files can be easily converted with programs like `samtools`


### Question 3
*What does the CIGAR string “26M2I78M” mean?*

- A CIGAR string describes the nature of bases in an alignment
- Corresponding directly to aligned bases read left to right, this CIGAR string reads "26 bases aligned, 2 bases inserted, 78 bases aligned" for a read of length $26+2+78=106$ bp
- Since M is an alignment match and could mean a sequence match or mismatch, its replacement with = for a sequence match and X for a sequence mismatch has been proposed
