# GSEA analysis pipeline
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(enrichplot)

# === Set working directory ===
setwd("/path/to/project/edgeR/")

# === Input sample comparisons ===
SMP <- c("ABISrtTA_v_HDFn", "rtTA_v_ABISrtTA", "rtTA_v_HDFn")
smp <- SMP[2]  # choose one comparison

# === Load data ===
data <- read.csv(paste0("prefix-", smp, "-all_genes.csv"))
head(data)

# === Prepare ranked gene list ===
geneList <- -data$logFC
names(geneList) <- data$Gene
geneList <- geneList[!duplicated(names(geneList))]
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[geneList != 0]

# === Get MSigDB GO sets ===
msig_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, gene_symbol) %>% dplyr::rename(ont = gs_name, gene = gene_symbol)
msig_cc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC") %>%
  dplyr::select(gs_name, gene_symbol) %>% dplyr::rename(ont = gs_name, gene = gene_symbol)
msig_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF") %>%
  dplyr::select(gs_name, gene_symbol) %>% dplyr::rename(ont = gs_name, gene = gene_symbol)

# === Run GSEA ===
gse_bp <- GSEA(geneList, TERM2GENE = msig_bp, pvalueCutoff = 1, pAdjustMethod = "none", maxGSSize = 2000)
gse_bp@result$go_sub <- "bp"

gse_cc <- GSEA(geneList, TERM2GENE = msig_cc, pvalueCutoff = 1, pAdjustMethod = "none", maxGSSize = 2000)
gse_cc@result$go_sub <- "cc"

gse_mf <- GSEA(geneList, TERM2GENE = msig_mf, pvalueCutoff = 1, pAdjustMethod = "none", maxGSSize = 2000)
gse_mf@result$go_sub <- "mf"

# === Save results ===
gse_tab <- rbind(gse_bp@result, gse_cc@result, gse_mf@result)
gse_tab$compare <- smp
write.csv(gse_tab, paste0(smp, "_go.csv"))

# === Extract top terms ===
top10_bp <- rbind(
  gse_bp@result[gse_bp@result$p.adjust < 0.05, ] %>% arrange(NES) %>% slice(1:5),
  gse_bp@result[gse_bp@result$p.adjust < 0.05, ] %>% arrange(desc(NES)) %>% slice(1:5)
)
top10_cc <- rbind(
  gse_cc@result[gse_cc@result$p.adjust < 0.05, ] %>% arrange(NES) %>% slice(1:5),
  gse_cc@result[gse_cc@result$p.adjust < 0.05, ] %>% arrange(desc(NES)) %>% slice(1:5)
)
top10_mf <- rbind(
  gse_mf@result[gse_mf@result$p.adjust < 0.05, ] %>% arrange(NES) %>% slice(1:5),
  gse_mf@result[gse_mf@result$p.adjust < 0.05, ] %>% arrange(desc(NES)) %>% slice(1:5)
)

top10_all <- rbind(top10_bp, top10_cc, top10_mf)
top10_all$group <- ifelse(top10_all$NES > 0, "Up", "Down")

# === Plot summary bar chart ===
p1 <- ggplot(top10_all, aes(NES, ID, fill = group)) +
  geom_col() +
  labs(x = "NES", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 9),
        legend.position = "right") +
  ggtitle(paste0(smp, " [p.adj<0.05]"))
ggsave(paste0(smp, "_topgo.png"), width = 10, height = 10, dpi = 300)

# === Venn diagram across comparisons ===
go_SUB <- c("bp", "mf", "cc")
for (i in go_SUB) {
  ABISrtTA_v_HDFn <- read.csv("ABISrtTA_v_HDFn_go.csv", row.names = 1)
  rtTA_v_ABISrtTA <- read.csv("rtTA_v_ABISrtTA_go.csv", row.names = 1)
  rtTA_v_HDFn <- read.csv("rtTA_v_HDFn_go.csv", row.names = 1)

  ABISrtTA_v_HDFn <- ABISrtTA_v_HDFn[ABISrtTA_v_HDFn$p.adjust < 0.05 & ABISrtTA_v_HDFn$go_sub == i, "Description"]
  rtTA_v_ABISrtTA <- rtTA_v_ABISrtTA[rtTA_v_ABISrtTA$p.adjust < 0.05 & rtTA_v_ABISrtTA$go_sub == i, "Description"]
  rtTA_v_HDFn <- rtTA_v_HDFn[rtTA_v_HDFn$p.adjust < 0.05 & rtTA_v_HDFn$go_sub == i, "Description"]

  vennList <- list(ABISrtTA_v_HDFn = ABISrtTA_v_HDFn,
                   rtTA_v_ABISrtTA = rtTA_v_ABISrtTA,
                   rtTA_v_HDFn = rtTA_v_HDFn)

  library(ggvenn)
  ggvenn(vennList,
         fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
         stroke_size = 0.5, set_name_size = 4) +
    ggtitle(i)

  ggsave(paste0(i, "_common_goterms.png"), width = 6, height = 6, dpi = 300)
}

# === Combine selected GO terms across comparisons ===
setwd("/path/to/project/edgeR/")
go_list <- read.csv("../go_list.csv", header = FALSE)
files <- list.files(pattern = "_go.csv")
smp <- gsub("_go.csv", "", files)

tocombine <- list()
for (i in 1:length(smp)) {
  tmp <- read.csv(files[i])
  tmp <- tmp[tmp$ID %in% go_list[, 1], ]
  tocombine[[i]] <- tmp
}
final_df <- Reduce(rbind, tocombine)
final_df$compare[final_df$compare == "rtTA_v_ABISrtTA"] <- "ABISrtTA_v_rtTA"

ggplot(final_df, aes(x = compare, y = Description, color = NES, size = -log10(pvalue))) +
  geom_point() +
  theme_bw() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  ggtitle("GO:BP")
ggsave("selected_go-bp.png", width = 11, height = 15, dpi = 300)

# === Subtraction analysis (remove baseline effect) ===
setwd("/path/to/project/")
go_list <- read.csv("go_list.csv", header = FALSE)

ABISrtTA <- read.csv("edgeR/ABISrtTA_v_HDFn_go.csv")
rtTA <- read.csv("edgeR/rtTA_v_HDFn_go.csv")

ABISrtTA_sub <- ABISrtTA[ABISrtTA$Description %in% go_list[, 1], ]
rtTA_sub <- rtTA[rtTA$Description %in% go_list[, 1], ]
common_go <- intersect(ABISrtTA_sub$Description, rtTA_sub$Description)

ABISrtTA_sub <- ABISrtTA_sub[ABISrtTA_sub$Description %in% common_go, c(3, 6)]
names(ABISrtTA_sub)[2] <- "ABISrtTA_NES"
rtTA_sub <- rtTA_sub[rtTA_sub$Description %in% common_go, c(3, 6)]
names(rtTA_sub)[2] <- "rtTA_NES"

merged_sub <- merge(ABISrtTA_sub, rtTA_sub, by = "Description")
merged_sub$substration <- merged_sub$ABISrtTA_NES - merged_sub$rtTA_NES
write.csv(merged_sub, "neuron-related_go_ABISrtTA_v_rtTA_after_remove_HDFn.csv", row.names = FALSE)

merged_sub$compare <- "ABISrtTA_v_rtTA_after_remove_HDFn"
merged_sub$group <- ifelse(merged_sub$substration > 0, "up", "down")

ggplot(merged_sub, aes(substration, Description, fill = group)) +
  geom_col() +
  labs(x = "NES", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 9),
        legend.position = "right") +
  ggtitle("Neuron-related GO terms: ABISrtTA_v_rtTA after removing HDFn")
ggsave("neuron-related_go_ABISrtTA_v_rtTA_after_remove_HDFn.png", width = 15, height = 10, dpi = 300)

