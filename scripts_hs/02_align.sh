#!/usr/bin/bash

#conda activate rna

# Move into your working directory
cd /path/to/project_directory
for i in `ls -1 trim/*1_paired.fq.gz | sed 's/_1_paired.fq.gz//' | sed 's/trim\///'`;
do
   mkdir align/${i}
   STAR --runThreadN 16 \
     --genomeDir /path/to/reference/star_index \
     --readFilesIn trim/${i}_1_paired.fq.gz trim/${i}_2_paired.fq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM \
     --outFileNamePrefix align/${i}/${i}
done
