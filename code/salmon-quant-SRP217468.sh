#!/bin/bash

for fn in data/datasets/SRP217468/SRR99049*/*.sra;
    do
        samp=`basename ${fn} .sra`
        dir=`dirname ${fn}`
        FILE="data/datasets/SRP217468/quants/${samp}_quant/quant.sf"
        echo "Checking file ${FILE}"

        if [ ! -f "$FILE" ]; then
            echo "Converting sample ${samp}"
            fastq-dump -I --split-files -O ${dir} ${fn}

            echo "Quantifying sample ${samp}"
            salmon quant -i ~/Work/genomes/hg19/salmon_v1.2.1_index -l A \
            -1 ${dir}/${samp}_1.fastq \
            -2 ${dir}/${samp}_2.fastq \
            -p 6 --validateMappings -o data/datasets/SRP217468/quants/${samp}_quant
        fi
done
