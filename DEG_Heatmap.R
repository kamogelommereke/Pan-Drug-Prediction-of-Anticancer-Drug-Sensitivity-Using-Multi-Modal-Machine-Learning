setwd(setwd("C:/DESKTOP/BIOINFORMATICS VALIDATION"))

###############################################################
## Differential Expression (Sensitive vs Resistant) — 12 Drugs
## End-to-end pipeline with explicit cleaning steps.
## Output: logFC matrix, FDR matrix, heatmap PDF
###############################################################

## =========================
## 0) Packages
## =========================
suppressPackageStartupMessages({
  library(tidyverse)     # readr, dplyr, tibble, ggplot2
  library(pheatmap)
})

## =========================
## 1) User parameters
## =========================
expr_file <- "Expression_clean_harmonized.csv"
auc_file  <- "GDSC_AUC_with_cellname_modelID.csv"
genes_file <- "TOP GENES.txt"

# 12 drugs of interest (names exactly as in GDSC table if possible)
drugs_12 <- c(
  "AZD8186", "Avagacestat", "BMS-536924", "Afatinib",
  "Foretinib", "AZD4547", "AZD7762", "AZD6738",
  "Gemcitabine", "Docetaxel", "Camptothecin", "5-Fluorouracil"
)

# Output prefixes
out_prefix <- "DE_12drugs"
dir.create("outputs", showWarnings = FALSE)

## =========================
## 2) Load & CLEAN the expression matrix
## =========================
message("Loading expression matrix: ", expr_file)
expr_raw <- readr::read_csv(expr_file, show_col_types = FALSE)

# ---- 2a) Choose an ID column (robustly)
id_candidates <- c("CellLineName_clean", "final_cellname", "CellLineName")
sample_col <- intersect(id_candidates, colnames(expr_raw))
if (length(sample_col) == 0) stop("No usable sample ID column found in expression file.")
expr <- expr_raw %>% dplyr::rename(Cell = !!sample_col[1])

# ---- 2b) Remove non-gene metadata columns (keep numeric only after 'Cell')
# Convert to data.frame to manage rownames easily
expr_df <- as.data.frame(expr)
# Clear any existing rownames (required before setting new ones)
rownames(expr_df) <- NULL

# (Optional) Inspect duplicates (print, but continue)
dups <- expr_df$Cell[duplicated(expr_df$Cell)]
if (length(dups) > 0) {
  message("Duplicate cell IDs detected (will be made unique): ",
          paste(unique(dups), collapse = ", "))
  expr_df$Cell <- make.unique(expr_df$Cell)   # e.g., TT -> TT, TT.1
}

# Set rownames from 'Cell' and drop the col
rownames(expr_df) <- expr_df$Cell
expr_df$Cell <- NULL

# Keep only numeric columns (genes)
is_num <- vapply(expr_df, is.numeric, logical(1L))
expr_df <- expr_df[, is_num, drop = FALSE]

# ---- 2c) Handle missing values & zero-variance genes (basic hygiene)
# Replace inf with NA, ensure numeric
expr_df[!is.finite(as.matrix(expr_df))] <- NA_real_
# Drop genes that are entirely NA
expr_df <- expr_df[, colSums(is.na(expr_df)) < nrow(expr_df), drop = FALSE]
# Optionally drop genes with near-zero variance across all lines (helps stability)
gene_var <- apply(expr_df, 2, stats::var, na.rm = TRUE)
expr_df <- expr_df[, gene_var > 0, drop = FALSE]

message("Expression matrix cleaned: ",
        nrow(expr_df), " samples × ", ncol(expr_df), " genes")

## =========================
## 3) Load & CLEAN your gene list (50 genes)
## =========================
message("Loading gene list: ", genes_file)
genes_raw <- readLines(genes_file, warn = FALSE)

# Remove ranks/numbers/commas; keep tokens that contain letters
genes_tokens <- strsplit(paste(genes_raw, collapse=" "), "\\s+")[[1]]
genes_tokens <- gsub(",", "", genes_tokens)
genes_clean <- genes_tokens[grepl("[A-Za-z]", genes_tokens)]

# Keep only genes present in expression
genes_found <- intersect(genes_clean, colnames(expr_df))
genes_missing <- setdiff(genes_clean, genes_found)
if (length(genes_found) == 0) stop("None of the requested genes found in expression table.")
if (length(genes_missing) > 0) {
  message("Missing genes (skipped): ", paste(genes_missing, collapse = ", "))
}

# Subset expression to your genes
expr_sub <- expr_df[, genes_found, drop = FALSE]

# Log1p transform (robust if counts/TPM-like; harmless if already logged)
expr_log <- log1p(expr_sub)

## =========================
## 4) Load & CLEAN the GDSC AUC table
## =========================
message("Loading AUC table: ", auc_file)
auc_raw <- readr::read_csv(auc_file, show_col_types = FALSE)

# ---- 4a) Choose a cell ID column that matches expression rownames
auc_id_candidates <- c("final_cellname", "CellLineName", "CellLineName_clean")
auc_id_col <- intersect(auc_id_candidates, colnames(auc_raw))
if (length(auc_id_col) == 0) stop("No usable cell ID column found in AUC file.")
auc <- auc_raw %>% dplyr::rename(Cell = !!auc_id_col[0 + 1]) %>%
  dplyr::mutate(Cell = as.character(Cell))

# ---- 4b) Keep only the 12 drugs of interest that actually exist in the AUC table
available_drugs <- intersect(drugs_12, colnames(auc))
if (length(available_drugs) == 0) stop("None of the 12 drugs were found in the AUC file.")
if (length(available_drugs) != length(drugs_12)) {
  message("These drugs were not found and will be skipped: ",
          paste(setdiff(drugs_12, available_drugs), collapse = ", "))
}
auc_sel <- auc %>% dplyr::select(Cell, dplyr::all_of(available_drugs))

# Coerce AUC columns to numeric and standardize Cell strings
auc_sel <- auc_sel %>%
  mutate(Cell = trimws(Cell)) %>%
  mutate(across(all_of(available_drugs), ~ suppressWarnings(as.numeric(.x))))

# Optional: drop rows with missing Cell
auc_sel <- auc_sel %>% dplyr::filter(!is.na(Cell) & Cell != "")

## =========================
## 5) Build Sensitive vs Resistant groups (per drug; tertiles of AUC)
## =========================
# Helper function: given a numeric vector, return bottom and top tertile thresholds
get_tertiles <- function(x) {
  if (all(is.na(x))) return(c(NA_real_, NA_real_))
  qs <- stats::quantile(x, probs = c(0.33, 0.67), na.rm = TRUE, names = FALSE)
  return(qs)
}

## =========================
## 6) Differential expression per drug
## =========================
genes <- colnames(expr_log)  # (genes_found)
logFC_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(available_drugs),
                    dimnames = list(genes, available_drugs))
FDR_mat   <- logFC_mat

group_sizes <- tibble(Drug = character(), n_sensitive = integer(), n_resistant = integer())

for (drug in available_drugs) {
  message("Processing ", drug, " ...")
  
  df <- auc_sel %>% dplyr::select(Cell, AUC = all_of(drug)) %>% dplyr::filter(!is.na(AUC))
  if (nrow(df) < 10) {  # safety: too few AUCs
    message("  Skipped: <", drug, "> has too few non-missing AUCs.")
    next
  }
  
  # Tertiles
  th <- get_tertiles(df$AUC)
  if (any(is.na(th))) {
    message("  Skipped: <", drug, "> tertiles could not be computed.")
    next
  }
  
  sens_ids <- df %>% dplyr::filter(AUC <= th[1]) %>% pull(Cell) %>% as.character()
  resi_ids <- df %>% dplyr::filter(AUC >= th[2]) %>% pull(Cell) %>% as.character()
  
  # Intersect with expression samples (row names)
  sens_ids <- intersect(sens_ids, rownames(expr_log))
  resi_ids <- intersect(resi_ids, rownames(expr_log))
  
  group_sizes <- bind_rows(
    group_sizes,
    tibble(Drug = drug, n_sensitive = length(sens_ids), n_resistant = length(resi_ids))
  )
  
  if (length(sens_ids) < 2 || length(resi_ids) < 2) {
    message("  Skipped: <", drug, "> not enough samples after alignment (S=", length(sens_ids),
            ", R=", length(resi_ids), ").")
    next
  }
  
  Xs <- expr_log[sens_ids, genes, drop = FALSE]
  Xr <- expr_log[resi_ids, genes, drop = FALSE]
  
  # Mean difference (Sensitive − Resistant)
  logFC <- colMeans(Xs, na.rm = TRUE) - colMeans(Xr, na.rm = TRUE)
  
  # Welch t-test per gene (vectorized via sapply)
  pvals <- sapply(genes, function(g) {
    xs <- Xs[, g]
    xr <- Xr[, g]
    # guard against all-NA
    if (all(is.na(xs)) || all(is.na(xr))) return(1)
    res <- tryCatch(t.test(xs, xr, var.equal = FALSE)$p.value, error = function(e) 1)
    if (!is.finite(res)) res <- 1
    res
  })
  
  # FDR per drug
  fdr <- p.adjust(pvals, method = "BH")
  
  logFC_mat[, drug] <- logFC
  FDR_mat[, drug]   <- fdr
}

# Save group sizes
readr::write_csv(group_sizes, file.path("outputs", paste0(out_prefix, "_group_sizes_per_drug.csv")))

## =========================
## 7) Save matrices
## =========================
logFC_df <- as.data.frame(logFC_mat)
FDR_df   <- as.data.frame(FDR_mat)

readr::write_csv(logFC_df, file.path("outputs", paste0(out_prefix, "_logFC_matrix.csv")))
readr::write_csv(FDR_df,   file.path("outputs", paste0(out_prefix, "_FDR_matrix.csv")))

message("Saved matrices to 'outputs/'")

## =========================
## 8) Heatmap of effect sizes (logFC)
## =========================
# Remove rows (genes) that are all NA
keep_genes <- rowSums(is.finite(as.matrix(logFC_mat))) > 0
logFC_hm <- logFC_mat[keep_genes, , drop = FALSE]

# Order genes optionally by mean |effect|
mean_abs <- apply(abs(logFC_hm), 1, mean, na.rm = TRUE)
logFC_hm <- logFC_hm[order(mean_abs, decreasing = TRUE), , drop = FALSE]

# Plot
pheatmap::pheatmap(
  logFC_hm,
  color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(255),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  legend = TRUE,
  main = "Differential Expression (Sensitive − Resistant) — 12 Drugs",
  fontsize_row = 7,
  fontsize_col = 9,
  border_color = NA,
  filename = file.path("outputs", paste0(out_prefix, "_logFC_heatmap.pdf")),
  width = 9, height = 12
)

message("Heatmap saved to 'outputs/", out_prefix, "_logFC_heatmap.pdf'")

## =========================
## 9) (Optional) Dot plot: color = logFC, size = −log10(FDR)
## =========================
# Uncomment this section if you want a quick dot plot summary

# library(ggplot2)
# df_long <- expand.grid(Gene = rownames(logFC_mat), Drug = colnames(logFC_mat), stringsAsFactors = FALSE) %>%
#   mutate(logFC = as.numeric(mapply(function(g, d) logFC_mat[g, d], Gene, Drug)),
#          FDR   = as.numeric(mapply(function(g, d) FDR_mat[g, d],   Gene, Drug))) %>%
#   filter(is.finite(logFC), is.finite(FDR))
#
# df_long$neglog10FDR <- -log10(pmax(df_long$FDR, 1e-300))
# L <- max(abs(df_long$logFC), na.rm = TRUE)
#
# gg <- ggplot(df_long, aes(Drug, Gene)) +
#   geom_point(aes(color = logFC, size = neglog10FDR), stroke = 0.2) +
#   scale_color_gradient2(low = "#3b4cc0", high = "#b40426", mid = "white", midpoint = 0, limits = c(-L, L)) +
#   scale_size(range = c(0.5, 4)) +
#   theme_minimal(base_size = 10) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "Drug–Gene DE: color = logFC (S−R), size = −log10(FDR)")
#
# ggsave(file.path("outputs", paste0(out_prefix, "_dotplot.pdf")), gg, width = 10, height = 12)

###############################################################
## End

dim(logFC_hm)
summary(as.vector(logFC_hm))
group_sizes
logFC_mat[1:5, ]
length(intersect(
  rownames(expr_log),
  auc_sel$Cell
))
## =========================
## 8) Heatmap of effect sizes (logFC) — PNG + PDF
## =========================

# Remove rows (genes) that are all NA
keep_genes <- rowSums(is.finite(as.matrix(logFC_mat))) > 0
logFC_hm <- logFC_mat[keep_genes, , drop = FALSE]

# Optional: Order genes by mean absolute effect size
mean_abs <- apply(abs(logFC_hm), 1, mean, na.rm = TRUE)
logFC_hm <- logFC_hm[order(mean_abs, decreasing = TRUE), , drop = FALSE]

# ---- SAVE AS PDF ----
pheatmap::pheatmap(
  logFC_hm,
  color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(255),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  legend = TRUE,
  main = "Differential Expression (Sensitive − Resistant) — 12 Drugs",
  fontsize_row = 7,
  fontsize_col = 9,
  border_color = NA,
  filename = file.path("outputs", paste0(out_prefix, "_logFC_heatmap.pdf")),
  width = 9,
  height = 12
)

# ---- SAVE AS PNG ----
png(
  filename = file.path("outputs", paste0(out_prefix, "_logFC_heatmap.png")),
  width = 1800, height = 2400, res = 300
)
pheatmap::pheatmap(
  logFC_hm,
  color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(255),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  legend = TRUE,
  main = "Differential Expression (Sensitive − Resistant) — 12 Drugs",
  fontsize_row = 7,
  fontsize_col = 9,
  border_color = NA
)
dev.off()

message("Saved heatmap as PDF and PNG in outputs/")
###############################################################

## =========================
## Heatmap of effect sizes (logFC) — PNG + PDF
## =========================

# Remove rows (genes) that are all NA
keep_genes <- rowSums(is.finite(as.matrix(logFC_mat))) > 0
logFC_hm <- logFC_mat[keep_genes, , drop = FALSE]

# Optional: Order genes by mean absolute effect size
mean_abs <- apply(abs(logFC_hm), 1, mean, na.rm = TRUE)
logFC_hm <- logFC_hm[order(mean_abs, decreasing = TRUE), , drop = FALSE]

# ---- SAVE AS PDF ----
pheatmap::pheatmap(
  logFC_hm,
  color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(255),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  legend = TRUE,
  main = "Differential Expression Heatmap",
  fontsize_row = 12,
  fontsize_col = 12,
  border_color = NA,
  filename = file.path("outputs", paste0(out_prefix, "_logFC_heatmap.pdf")),
  width = 9,
  height = 12
)

# ---- SAVE AS PNG ----
png(
  filename = file.path("outputs", paste0(out_prefix, "_logFC_heatmap.png")),
  width = 1800, height = 2400, res = 300
)
pheatmap::pheatmap(
  logFC_hm,
  color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(255),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  legend = TRUE,
  main = "Differential Expression Heatmap",
  fontsize_row = 7,
  fontsize_col = 9,
  border_color = NA
)
dev.off()

message("Saved PDF + PNG heatmaps in outputs/")

getwd()
list.files("outputs")