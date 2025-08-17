#!/usr/bin/bash

# Set Java temporary directory
export _JAVA_OPTIONS="-Djava.io.tmpdir=/path/to/tmpdir"

# Move into project directory
cd /path/to/project_directory
mkdir -p trim

# Loop through raw fastq files
for i in `ls -1 raw/*1.fq.gz | sed 's/_1.fq.gz//' | sed 's/raw\///'`;
do
    trimmomatic PE -threads 16 \
        raw/${i}_1.fq.gz raw/${i}_2.fq.gz \
        trim/${i}_1_paired.fq.gz trim/${i}_1_unpaired.fq.gz \
        trim/${i}_2_paired.fq.gz trim/${i}_2_unpaired.fq.gz \
        ILLUMINACLIP:/path/to/adapters/TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
done

# Remove unpaired reads
rm trim/*_unpaired.fq.gz
