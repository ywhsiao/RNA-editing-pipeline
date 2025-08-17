# hs RNAseq using TPM and pairwise T test
rm(list = ls())
library(dplyr)
library(scran)

# === Define batch name and working directory ===
batch_name <- "batch1"
setwd(paste0("/path/to/project/", batch_name, "_hs/quan/quan_gene/"))

# === Metadata ===
files <- list.files(pattern = ".genes.result")
smp <- gsub("-LF.*", "", files)
condition <- gsub("\\-[0-9]", "", smp)  # adjust depending on naming convention
metadata <- cbind(as.data.frame(smp), as.data.frame(condition))
names(metadata) <- c("sample_id", "treatment")
metadata

# === TPM expression matrix ===
exp_list <- list()
for (i in 1:length(files)) {
  exp_list[[i]] <- read.csv(files[i], sep = "\t")[, c(1, 6)]
}
exp <- as.data.frame(Reduce(cbind, exp_list))
row.names(exp) <- exp[, 1]
exp_tpm <- exp[, -seq(1, ncol(exp), 2)]
names(exp_tpm) <- smp
exp_tpm$gene <- gsub("\\..*", "", row.names(exp_tpm))

# === Annotation (ENSEMBL → SYMBOL) ===
require(org.Hs.eg.db)
mapping <- mapIds(org.Hs.eg.db,
                  keys = exp_tpm$gene,
                  column = "SYMBOL",
                  keytype = "ENSEMBL")
final_mapping <- data.frame(SYMBOL = mapping, ENSEMBL = names(mapping))
exp_tpm_annot <- merge(final_mapping, exp_tpm, by.x = "ENSEMBL", by.y = "gene")
exp_tpm_annot <- exp_tpm_annot[, -1]
names(exp_tpm_annot)[1] <- "gene_symbol"

exp_annot_final <- exp_tpm_annot %>%
  group_by(gene_symbol) %>%
  summarise_all("mean") %>% as.data.frame()
exp_annot_final <- exp_annot_final[!is.na(exp_annot_final$gene_symbol), ]
row.names(exp_annot_final) <- exp_annot_final$gene_symbol
exp_annot_final <- exp_annot_final[, -1]
exp_annot_final <- exp_annot_final[rowSums(exp_annot_final[]) > 1, ]

# Save processed TPM matrices
write.csv(exp_annot_final, paste0(batch_name, "_tpm_notlog2.csv"))
exp_tpm_annot_final <- log2(exp_annot_final + 1)
write.csv(exp_tpm_annot_final, paste0(batch_name, "_tpm.csv"))

# === Boxplot ===
library(RColorBrewer)
png("exp_tpm_boxplot.png", width = 2500, height = 1500, res = 300)
boxplot(exp_tpm_annot_final, las = 2, main = "")
title(ylab = "log2(TPM+1)")
dev.off()

# === PCA ===
library(ggfortify)
df_final_t <- t(exp_tpm_annot_final)
df_meta <- cbind(metadata, df_final_t)
pca_res <- prcomp(df_meta[, -(1:2)])
autoplot(pca_res, label = TRUE, data = df_meta, colour = "treatment", label.size = 5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("exp_tpm_pca.png", width = 8, height = 8, dpi = 300)

# === Pairwise t-test ===
treatment <- metadata$treatment
sce <- SingleCellExperiment(list(counts = exp_tpm_annot_final),
                            colData = DataFrame(label = treatment))
out <- pairwiseTTests(counts(sce), groups = treatment)
saveRDS(out, "pairwise_de.rds")

# === Extract DE results per batch ===
work_dir <- paste0("/path/to/project/", batch_name, "_hs/quan/quan_gene/")
dat <- readRDS(paste0(work_dir, "pairwise_de.rds"))
for (k in c(3, 5)) {
  de <- as.data.frame(dat$statistics[[k]])
  write.csv(de, paste0(work_dir, dat$pairs[k, 1], "_vs_", dat$pairs[k, 2], "_de.csv"))
}

# === Volcano plots ===
library(ggplot2)
library(ggrepel)
data_lst <- list()
for (k in c(3, 5)) {
  data <- as.data.frame(out$statistics[[k]])
  data$logp <- -log10(data$p.value)
  data_lst[[k]] <- data
}

k <- 5
ggplot(data_lst[[k]], aes(logFC, logp)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = logp, color = logp)) +
  scale_color_gradientn(colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw() +
  guides(col = guide_colourbar(title = "-Log10_p-value"),
         size = "none") +
  xlab("Log2FC") +
  ylab("-Log10(p-value)") +
  ggtitle(paste0(out$pairs[k, 1], "_vs_", out$pairs[k, 2]))
ggsave(paste0(out$pairs[k, 1], "_vs_", out$pairs[k, 2], "_volcano.pdf"), height = 6, width = 8, dpi = 1200)

# === GO Enrichment ===
library(msigdbr)
library(clusterProfiler)
dir.create(paste0("/path/to/project/", batch_name, "_hs/quan/quan_gene/go"))
for (k in c(3, 5)) {
  data <- as.data.frame(out$statistics[[k]])
  geneList <- data$logFC
  names(geneList) <- gsub("\\|.*", "", row.names(data))
  geneList <- geneList[!duplicated(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)
  msig_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
    dplyr::select(gs_name, gene_symbol)
  gse_bp <- GSEA(geneList, TERM2GENE = msig_bp, pvalueCutoff = 1)
  write.csv(gse_bp@result,
            paste0("/path/to/project/", batch_name, "_hs/quan/quan_gene/go/",
                   out$pairs[k, 1], "_vs_", out$pairs[k, 2], "_allgo.csv"))
}

# === Venn Diagram ===
library(eulerr)
deg_lst <- list()
lst_name <- c()
for (k in c(3, 5)) {
  data <- as.data.frame(out$statistics[[k]])
  deg_lst[[k]] <- row.names(data[abs(data$logFC) >= 1 & data$p.value < 0.05, ])
  lst_name[k] <- paste0(out$pairs[k, 1], "_vs_", out$pairs[k, 2])
}
names(deg_lst) <- lst_name
fit1 <- euler(deg_lst)
pdf("area-proportional_venn.pdf")
plot(fit1, legend = TRUE,
     quantities = list(type = c("counts", "percent"), font = 3, round = 2, cex = 0.8))
dev

