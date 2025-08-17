#!/usr/bin/bash
  
#conda activate base

cd /mnt/053768d5-99c6-47ed-b61a-7c1f94831986/yiwen/mut_call/batch3

for i in `ls -1 trim/*1_paired.fq.gz | sed 's/\_1_paired.fq.gz//'| sed 's/trim\///'`;
do
	samtools index -@ 24 align/${i}/${i}Aligned.sortedByCoord.out.bam
done
					    
