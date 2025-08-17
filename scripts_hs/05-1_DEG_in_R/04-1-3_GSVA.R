# GSVA for sample-based study
# Set your project directory
setwd("/path/to/project/")

# ---- Load & transform expression ----
exp <- read.csv("gene_counts.txt", sep = "\t", row.names = 1)
cpm <- t(t(exp) / colSums(exp) * 1e6)
cpm <- log2(cpm + 1)

boxplot(cpm, ylab = expression('Log'[2]~'(CPM+1)'), las = 2, main = "CPM")

# ---- GSVA (GO: BP/CC/MF) ----
library(clusterProfiler)
library(msigdbr)
library(GSVA)

subcat <- c("BP", "CC", "MF")
for (subcategory in subcat) {
  GO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = subcategory)
  GO_df <- dplyr::select(GO, gs_name, gene_symbol, gs_exact_source, gs_subcat)
  go_list <- split(GO_df$gene_symbol, GO_df$gs_name)

  gsva_mat <- gsva(
    expr = as.matrix(cpm),
    gset.idx.list = go_list,
    kcdf = "Gaussian",           # "Gaussian" for logCPM/logRPKM/logTPM; "Poisson" for counts
    verbose = TRUE,
    parallel.sz = parallel::detectCores()
  )

  write.csv(as.data.frame(gsva_mat), paste0("GSVA_GO-", subcategory, ".csv"))
}

# ---- Hallmark (MSigDB H) ----
H <- msigdbr(species = "Homo sapiens", category = "H")
H_df <- dplyr::select(H, gs_name, gene_symbol, gs_exact_source, gs_subcat)
H_list <- split(H_df$gene_symbol, H_df$gs_name)

gsva_mat <- gsva(
  expr = as.matrix(cpm),
  gset.idx.list = H_list,
  kcdf = "Gaussian",
  verbose = TRUE,
  parallel.sz = parallel::detectCores()
)
write.csv(as.data.frame(gsva_mat), "GSVA_Hallmark.csv")

# ---- Reactome (C2:REACTOME) ----
R <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
R_df <- dplyr::select(R, gs_name, gene_symbol, gs_exact_source, gs_subcat)
R_list <- split(R_df$gene_symbol, R_df$gs_name)

gsva_mat <- gsva(
  expr = as.matrix(cpm),
  gset.idx.list = R_list,
  kcdf = "Gaussian",
  verbose = TRUE,
  parallel.sz = parallel::detectCores()
)
write.csv(as.data.frame(gsva_mat), "GSVA_Reactome.csv")

# ---- KEGG (C2:KEGG) ----
K <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
K_df <- dplyr::select(K, gs_name, gene_symbol, gs_exact_source, gs_subcat)
K_list <- split(K_df$gene_symbol, K_df$gs_name)

gsva_mat <- gsva(
  expr = as.matrix(cpm),
  gset.idx.list = K_list,
  kcdf = "Gaussian",
  verbose = TRUE,
  parallel.sz = parallel::detectCores()
)
write.csv(as.data.frame(gsva_mat), "GSVA_KEGG.csv")

# ---- Example heatmap for selected GO:BP terms ----
dat <- read.csv("GSVA_GO-BP.csv", row.names = 1)
subterms <- c(
  "GOBP_EXCITATORY_CHEMICAL_SYNAPTIC_TRANSMISSION",
  "GOBP_MAINTENANCE_OF_POSTSYNAPTIC_SPECIALIZATION_STRUCTURE",
  "GOBP_NEUROMUSCULAR_SYNAPTIC_TRANSMISSION",
  "GOBP_NEUROTRANSMITTER_LOADING_INTO_SYNAPTIC_VESICLE",
  "GOBP_REGULATION_OF_SHORT_TERM_NEURONAL_SYNAPTIC_PLASTICITY"
)
dat <- dat[row.names(dat) %in% subterms, ]

pheatmap::pheatmap(
  dat,
  fontsize_row = 8,
  height = 10,
  width = 18,
  show_colnames = TRUE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  filename = "GSVA_go_heatmap.pdf"
)

