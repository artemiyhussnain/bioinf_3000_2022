#!/bin/bash

echo Running NGS alignment pipeline
echo Expected running time ~1h
echo -e '\nDownloading and indexing reference genome'
bash script1.sh
echo -e '\nDownloading, trimming, and aligning NGS sequencing data'
bash script2.sh
echo -e '\nGenerating alignment summary'
bash script3.sh
echo -e '\nAll Done'