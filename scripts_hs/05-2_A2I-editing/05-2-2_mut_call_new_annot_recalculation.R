# ================================
# RNA Editing Analysis Pipeline
# Version: 2024.04
# ================================

# ---------- (1) Global mutation patterns ----------
setwd("<PROJECT_DIR>/batchX/annot2")   # <-- replace with your batch directory
files <- list.files(pattern = '.txt')
sample_names <- gsub('\\-LF.*','',files)
tissue <- gsub('\\-[0-9]','',sample_names)

for (i in 1:length(files)){
  df <- read.csv(files[i], sep = '\t')
  df_sub1 <- df[df$Ref!='-',]
  df_sub2 <- df_sub1[df_sub1$Alt!='-',]
  
  df_final <- df_sub2[,c(1,2,4,5,45)]
  df_sub <- df_final[,-5]
  names(df_sub) <- c('CHROM', 'POS', 'REF', 'ALT')
  df_sub$ID <- '.'
  df_sub$QUAL <- '.'
  df_sub$FILTER <- 'PASS'
  df_sub$INFO <- paste0('DP=',df_final$Otherinfo3)
  df_sub$FORMAT <- 'GT'
  df_sub$smp <- '0/1'
  smp <- sample_names[i]
  names(df_sub)[10] <- smp
  df_sub <- df_sub[,c(1,2,5,3,4,6,7,8,9,10)]
  
  CON <- file(paste0(smp,'_clean.txt'), "w")
  writeLines('##fileformat=VCFv4.2', CON)
  writeLines('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', CON)
  writeLines(paste0('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t',smp), CON)
  write.table(df_sub, CON, append=T, quote=FALSE, sep="\t", col.names = F, row.names = F)
  close(CON)
  file.rename(paste0(smp,'_clean.txt'), paste0(smp,'_clean.vcf'))
}

rm(list = ls())
library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

vcf_files <- list.files(pattern = "_clean.vcf$")
sample_names <- gsub('_clean.vcf','',vcf_files)
tissue <- gsub('\\-[0-9]','',sample_na_

