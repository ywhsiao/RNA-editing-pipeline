# RNA-editing-pipeline

This pipeline integrates bulk RNA-seq preprocessing, alignment, expression quantification, differential expression analysis, A-to-I RNA-editing profiling, and functional annotation.

**Highlight of this pipeline**

This pipeline not only quantifies transcriptional and RNA-editing changes, but also interprets them in a biological context by mapping results onto pathways and gene sets (via GSEA/GSVA).

![Preview](https://raw.githubusercontent.com/ywhsiao/RNA-editing-pipeline/master/png/potential_data_vis.png)
## Figure 1. Potential downstream analysis

1. Preprocessing and Alignment

- Input: Paired-end FASTQ files.
  
- Quality control and trimming were applied to remove adapters and low-quality reads.
  
- Alignment: Cleaned reads were aligned to the Homo sapiens reference genome (hg38) using STAR, generating:

- Genome-aligned BAMs (Aligned.sortedByCoord.out.bam)

- Transcriptome-aligned BAMs (Aligned.toTranscriptome.out.bam)

2. Gene Expression Quantification

- RSEM was applied to transcriptome-aligned BAMs.

- Output includes raw counts, TPMs, and FPKMs.

3. Differential Expression Analysis (DEGs)

- Count matrices were imported into customized R script for normalization and statistical testing.

- Differentially expressed genes were identified across conditions, controlling for multiple testing.

4. Variant Calling for RNA-editing Candidates

- Genome-aligned BAMs were preprocessed with samtools (read group addition, indexing).

- Mutect2 (GATK) was used for somatic SNV calling, with germline resource and panel of normals.

- Candidate SNVs were filtered and annotated.

5. Variant Annotation

- Annotation performed using ANNOVAR, incorporating databases such as:

  - refGene for gene context

  - dbSNP (avsnp150) for common SNPs

  - ClinVar for pathogenic variants

  - REVEL/functional predictors

- Output: multi-annotated variant tables (*_multianno.txt).

6. RNA-editing Site Identification

- Filtering pipeline:

  - Remove indels, retain only SNVs.

  - Exclude variants with missing REF/ALT bases.

  - Restrict to A/T/G/C bases.

  - Highlight A-to-I conversions (A>G, T>C).

- Cleaned VCF files were generated for downstream analysis.

7. Mutational Signature & A-to-I Profiling

- VCFs were imported into MutationalPatterns in R.

- 96-class mutational signatures were generated and relabeled to highlight RNA-editing–specific substitutions.

- Context-specific editing distributions were computed (5′/3′ neighbors, mutation categories).

8. Functional Annotation (GSEA & GSVA)

- Gene Set Enrichment Analysis (GSEA):

  - DEGs (ranked by log2 fold-change) were analyzed against MSigDB collections (e.g., Hallmark, KEGG, Reactome).

  - Identifies pathways enriched among up- and down-regulated genes.

- Gene Set Variation Analysis (GSVA):

  - TPM expression matrix was used for sample-level pathway enrichment.

  - GSVA scores reflect relative activation of pathways across tissues/conditions.

- Allows integration of expression and editing profiles to explore whether RNA-editing–associated genes cluster within specific pathways.

 
  
