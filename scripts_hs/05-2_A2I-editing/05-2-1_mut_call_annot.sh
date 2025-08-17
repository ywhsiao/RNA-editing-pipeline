#!/usr/bin/bash

# Move into project directory
cd /path/to/project_directory/
mkdir -p annot2
mkdir -p vcf2

for i in `ls -1 align/*Aligned.sortedByCoord.out.bam | sed 's/Aligned.sortedByCoord.out.bam//' | sed 's/align\///'`;
do
     # Add read groups & index
     samtools addreplacerg -@16 -r '@RG\tID:samplename\tSM:samplename' \
         align/${i}Aligned.sortedByCoord.out.bam \
         -o align/${i}Aligned.sortedByCoord.out.revised.bam
     samtools index -@16 align/${i}Aligned.sortedByCoord.out.revised.bam

     # Run Mutect2
     java -Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true \
          -jar /path/to/tools/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar Mutect2 \
          -R /path/to/reference/hg38_all.fa \
          -I align/${i}Aligned.sortedByCoord.out.revised.bam \
          --germline-resource /path/to/resources/af-only-gnomad.hg38.vcf.gz \
          --panel-of-normals /path/to/resources/1000g_pon.hg38.vcf.gz \
          --disable-read-filter MappingQualityAvailableReadFilter \
          -O vcf2/${i}_single_sample.vcf.gz

     # Extract SNPs
     gunzip -k vcf2/${i}_single_sample.vcf.gz
     bcftools view -v snps vcf2/${i}_single_sample.vcf > vcf2/${i}_single_sample_SNP.vcf

     # Annotate with ANNOVAR
     perl /path/to/tools/annovar/table_annovar.pl \
          vcf2/${i}_single_sample_SNP.vcf \
          /path/to/tools/annovar/humandb \
          -buildver hg38 \
          -out ./ \
          -remove \
          -protocol refGene,avsnp150,ljb26_all,clinvar_20221231,revel \
          -operation g,f,f,f,f \
          -nastring . \
          -vcfinput \
          --outfile annot2/${i}_SNP_anno
done
