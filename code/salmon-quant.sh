#!/bin/bash
for fn in data/SRP125125/SRR6298*;
do
samp=`basename ${fn}`
echo "Converting sample ${samp}"
fastq-dump -I --split-files -O ${fn} ${fn}/${samp}.sra

echo "Quantifying sample ${samp}"
salmon quant -i ~/Work/genomes/hg19/salmon_v1.2.1_index -l A \
          -1 ${fn}/${samp}_1.fastq \
          -2 ${fn}/${samp}_2.fastq \
          -p 4 --validateMappings -o data/SRP125125/quants/${samp}_quant
done
