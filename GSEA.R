###############################################################
## GSEA: Hallmark & Cancer Hallmark Pathways
## Uses msigdbr + clusterProfiler
###############################################################

library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(fgsea)

## -----------------------------------------
## 1) Load DE results (logFC & FDR)
## -----------------------------------------
logFC_mat <- read.csv("outputs/DE_12drugs_logFC_matrix.csv",
                      row.names = 1,
                      check.names = FALSE)

logFC_mat <- as.matrix(logFC_mat)

list.files()
list.files("outputs")

# Load FDR matrix WITHOUT rownames
FDR_mat <- read.csv("outputs/DE_12drugs_FDR_matrix.csv",
                    check.names = FALSE,
                    stringsAsFactors = FALSE)

# Identify the gene-name column (first column)
gene_col <- colnames(FDR_mat)[1]

# ---- FIX: Force gene name column to character ----
FDR_mat[[gene_col]] <- as.character(FDR_mat[[gene_col]])

# Make gene names unique
FDR_mat[[gene_col]] <- make.unique(FDR_mat[[gene_col]])

# Move to rownames
rownames(FDR_mat) <- FDR_mat[[gene_col]]

# Drop gene-name column
FDR_mat[[gene_col]] <- NULL

# Convert to matrix
FDR_mat <- as.matrix(FDR_mat)

## -----------------------------------------
## 2) Choose pathways from MSigDB
## -----------------------------------------

# Hallmark gene sets (H)
msig_h <- msigdbr(species = "Homo sapiens", category = "H")

# Cancer-related collections:
# C6 = Oncogenic Signatures (pre-2020)
# C7 = Immunologic Signatures, sometimes used in tumor biology
msig_cancer <- msigdbr(species = "Homo sapiens", category = "C6")

# Prepare GMT-like list for fgsea
hallmark_list <- msig_h %>% 
  split(x = .$gene_symbol, f = .$gs_name)

cancerhallmark_list <- msig_cancer %>% 
  split(x = .$gene_symbol, f = .$gs_name)

## -----------------------------------------
## 3) Function to create ranked gene list
## -----------------------------------------
rank_genes <- function(logfc_vec, fdr_vec) {
  
  # 1) Replace NA with 0
  logfc_vec[is.na(logfc_vec)] <- 0
  fdr_vec[is.na(fdr_vec)] <- 1
  
  # 2) Build composite ranking metric
  # Signed score = logFC * -log10(FDR)
  score <- logfc_vec * (-log10(fdr_vec + 1e-300))
  
  # 3) Sort decreasing
  score <- sort(score, decreasing = TRUE)
  
  return(score)
}

## -----------------------------------------
## 4) Run GSEA for each drug
## -----------------------------------------

drugs <- colnames(logFC_mat)

# Create output directory
dir.create("GSEA_outputs", showWarnings = FALSE)

for (drug in drugs) {
  
  message("Running GSEA for: ", drug)
  
  # Extract gene vectors
  logfc_vec <- logFC_mat[, drug]
  fdr_vec   <- FDR_mat[, drug]
  
  ranks <- rank_genes(logfc_vec, fdr_vec)
  
  ## ---- Run FGSEA for Hallmark ----
  fg_h <- fgsea(
    pathways = hallmark_list,
    stats = ranks,
    nperm = 10000
  )
  
  ## ---- Run FGSEA for Cancer Hallmark ----
  fg_cancer <- fgsea(
    pathways = cancerhallmark_list,
    stats = ranks,
    nperm = 10000
  )
  
  # Save results
  write_csv(fg_h,     paste0("GSEA_outputs/", drug, "_Hallmark_GSEA.csv"))
  write_csv(fg_cancer, paste0("GSEA_outputs/", drug, "_CancerHallmark_GSEA.csv"))
  
  # Save top plots
  pdf(paste0("GSEA_outputs/", drug, "_Hallmark_GSEA_topplots.pdf"), width=8, height=10)
  plotEnrichment(hallmark_list[[fg_h$pathway[1]]], ranks)
  dev.off()
  
  pdf(paste0("GSEA_outputs/", drug, "_CancerHallmark_GSEA_topplots.pdf"), width=8, height=10)
  plotEnrichment(cancerhallmark_list[[fg_cancer$pathway[1]]], ranks)
  dev.off()
}

message("GSEA completed for all drugs!")