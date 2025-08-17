#!/usr/bin/bash
  
# Activate conda environment
conda activate rna

# Move into project directory
cd /path/to/project_directory
mkdir -p quan

# Loop through trimmed fastq files
for i in `ls -1 trim/*1_paired.fq.gz | sed 's/_1_paired.fq.gz//' | sed 's/trim\///'`;
do
    mkdir -p quan/${i}
    rsem-calculate-expression \
        --bam --estimate-rspd --no-bam-output \
        -p 16 --paired-end \
        align/${i}/${i}Aligned.toTranscriptome.out.bam \
        /path/to/reference/rsem_index/grch38 \
        quan/${i}/${i}
done
