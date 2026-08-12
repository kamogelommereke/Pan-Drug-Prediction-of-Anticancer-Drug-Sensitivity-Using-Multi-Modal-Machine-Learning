setwd(setwd("C:/DESKTOP/BIOINFORMATICS VALIDATION"))

###############################################################################
# FULL SCRIPT: 12 Drug Volcano Plots (Stacked on One Page)
# - Sensitive vs Resistant (bottom 25% vs top 25% AUC)
# - Blue/Green colors
# - Genes labeled (TOP GENES or fallback)
# - One combined PDF with all 12 plots stacked vertically
###############################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
`
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)     # for stacking plots
})

set.seed(1234)

###############################################################################
# INPUT FILES
###############################################################################
exp_path   <- "Expression_clean_harmonized.csv"
gdsc_path  <- "GDSC_AUC_with_cellname_modelID.csv"
model_path <- "Model.csv"
top_genes_file <- "TOP GENES.txt"

###############################################################################
# PARAMETERS
###############################################################################
auc_lo <- 0.25
auc_hi <- 0.75
alpha_fdr <- 0.05
logFC_cut <- 0.58    # ≈1.5 fold
min_group <- 5
label_top_n <- 15
png_res <- 300

###############################################################################
# EXACT DRUG LIST (from your GDSC file)
###############################################################################
focus_drugs <- c(
  "AZD8186", "Avagacestat", "BMS-536924", "Afatinib", "Foretinib",
  "AZD4547", "AZD7762", "AZD6738", "Gemcitabine",
  "Docetaxel", "Camptothecin", "5-Fluorouracil"
)

###############################################################################
# NORMALIZATION function
###############################################################################
normalize_name <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]", "", x)
  toupper(x)
}

first_non_na <- function(...) {
  vals <- c(...)
  vals[which(!is.na(vals))[1]]
}

###############################################################################
# LOAD TOP GENES
###############################################################################
top_raw <- readLines(top_genes_file, warn = FALSE)
extract_symbol <- function(line){
  toks <- unlist(strsplit(trimws(line), "\\s+"))
  toks[length(toks)]
}
TOP_GENES <- unique(na.omit(sapply(top_raw, extract_symbol)))

###############################################################################
# LOAD MODEL MAP
###############################################################################
model_df <- read_csv(model_path, show_col_types = FALSE) %>%
  mutate(
    StrippedCellLineName = if("StrippedCellLineName"%in%names(.)) StrippedCellLineName else NA,
    CellLineName         = if("CellLineName"%in%names(.)) CellLineName         else NA,
    NormalizedName       = coalesce(normalize_name(StrippedCellLineName),
                                    normalize_name(CellLineName))
  ) %>%
  select(ModelID, NormalizedName)

###############################################################################
# LOAD GDSC
###############################################################################
gdsc <- read_csv(gdsc_path, show_col_types = FALSE)

gdsc <- gdsc %>%
  mutate(
    cellname_modelID_clean = if("cellname_modelID"%in%names(.)) normalize_name(cellname_modelID) else NA,
    final_cellname_clean   = if("final_cellname"%in%names(.))   normalize_name(final_cellname)   else NA,
    CellLineName_clean     = if("CellLineName"%in%names(.))     normalize_name(CellLineName)     else NA
  ) %>%
  left_join(model_df, by="ModelID") %>%
  mutate(HarmonizedName = pmap_chr(
    list(cellname_modelID_clean, final_cellname_clean,
         CellLineName_clean, NormalizedName),
    ~ first_non_na(..1, ..2, ..3, ..4)
  ))

###############################################################################
# LOAD EXPRESSION
###############################################################################
exp <- read_csv(exp_path, show_col_types = FALSE)
exp <- exp %>%
  mutate(CellLineName_clean = normalize_name(CellLineName)) %>%
  distinct(CellLineName_clean, .keep_all = TRUE)

###############################################################################
# JOIN EXPRESSION + GDSC
###############################################################################
joined <- exp %>%
  inner_join(
    gdsc %>% select(HarmonizedName, all_of(focus_drugs)),
    by = c("CellLineName_clean" = "HarmonizedName")
  )

expr_cols <- setdiff(names(exp), c("CellLineName","CellLineName_clean"))

###############################################################################
# STORAGE FOR 12 VOLCANO PLOTS
###############################################################################
volcano_list <- list()

###############################################################################
# LOOP OVER ALL 12 DRUGS
###############################################################################
for (drug in focus_drugs) {
  
  message("Processing drug: ", drug)
  
  # 1) Compute groups
  auc <- joined[[drug]]
  qlo <- quantile(auc, auc_lo, na.rm = TRUE)
  qhi <- quantile(auc, auc_hi, na.rm = TRUE)
  
  joined2 <- joined %>%
    mutate(Group = ifelse(auc <= qlo, "Sensitive",
                          ifelse(auc >= qhi, "Resistant", NA))) %>%
    filter(!is.na(Group))
  
  if (length(unique(joined2$Group)) < 2 || min(table(joined2$Group)) < min_group) {
    message("  Skipped: insufficient samples.")
    next
  }
  
  joined2$Group <- factor(joined2$Group, levels = c("Sensitive","Resistant"))
  
  # 2) Expression matrix
  E <- t(as.matrix(joined2[, expr_cols]))
  
  # Remove zero variance genes
  E <- E[apply(E,1,var,na.rm=TRUE) > 0, , drop=FALSE]
  if (nrow(E) == 0) {
    message("  Skipped: zero variance genes.")
    next
  }
  
  # 3) LIMMA
  design <- model.matrix(~ Group, data = joined2)
  fit <- lmFit(E, design)
  fit <- eBayes(fit, trend = TRUE)
  
  coef_name <- "GroupResistant"
  if (!(coef_name %in% colnames(fit$coefficients))) {
    message("  Skipped: missing contrast.")
    next
  }
  
  tt <- topTable(fit, coef=coef_name, number=Inf, p.value=1, adjust.method="BH")
  if (nrow(tt)==0) {
    message("  Skipped: empty DE table.")
    next
  }
  
  tt$gene <- rownames(tt)
  tt <- tt %>%
    rename(log2FC = logFC, PValue=P.Value, FDR=adj.P.Val) %>%
    relocate(gene)
  
  tt$neglog10FDR <- -log10(pmax(tt$FDR,1e-300))
  tt$sig <- (tt$FDR < alpha_fdr & abs(tt$log2FC) >= logFC_cut)
  
  # 4) Label genes (TOP GENES OR fallback top N)
  label_df <- tt %>% filter(sig & gene %in% TOP_GENES)
  if (nrow(label_df)==0) {
    label_df <- tt %>%
      filter(sig) %>%
      slice_max(order_by = abs(log2FC), n = label_top_n, with_ties = FALSE)
  }
  
  # 5) Build volcano plot
  p <- ggplot(tt, aes(x=log2FC, y=neglog10FDR)) +
    geom_point(aes(color=sig), alpha=0.7, size=1.3) +
    scale_color_manual(values=c(`TRUE`="#1f77b4", `FALSE`="#2ca02c")) +
    geom_vline(xintercept=c(-logFC_cut, logFC_cut),
               linetype=2, color="grey55") +
    geom_hline(yintercept=-log10(alpha_fdr),
               linetype=2, color="grey55") +
    ggrepel::geom_text_repel(
      data=label_df, aes(label=gene), size=2.7,
      max.overlaps=Inf, min.segment.length=0
    ) +
    labs(
      title = drug,
      x = "log2FC (Resistant - Sensitive)",
      y = expression(-log10(FDR))
    ) +
    theme_minimal(base_size = 11)
  
  volcano_list[[drug]] <- p
}

###############################################################################
# STACK ALL 12 VOLCANO PLOTS INTO 1 PAGE (vertical)
###############################################################################
combined <- cowplot::plot_grid(
  plotlist = volcano_list,
  ncol = 1,
  align = "v",
  rel_heights = rep(1, length(volcano_list))
)

# Save combined figure
ggsave("ALL_12_DRUGS_STACKED_VOLCANO.pdf",
       combined, width = 9, height = 45)  # tall PDF

message("=============================================")
message(" ALL 12 DRUG VOLCANO PLOTS SAVED IN ONE PDF ")
message(" File: ALL_12_DRUGS_STACKED_VOLCANO.pdf      ")
message("=============================================")




###############################################################################
# 12-Drug Volcano Grid (4 per row)
# - Sensitive vs Resistant (bottom 25% vs top 25% AUC)
# - One overall title (no per-plot titles)
# - Gene labels (TOP GENES; fallback to top-|log2FC| significant)
# - BLUE = significant; GREEN = non-significant
# - Saves PDF + PNG of the page
# - Exports one CSV per drug + one combined CSV with all plotted points
###############################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)   # for grid + overall title
})

set.seed(1234)

# -------------------------
# INPUT FILES
# -------------------------
exp_path        <- "Expression_clean_harmonized.csv"
gdsc_path       <- "GDSC_AUC_with_cellname_modelID.csv"
model_path      <- "Model.csv"
top_genes_file  <- "TOP GENES.txt"

# -------------------------
# PARAMETERS
# -------------------------
auc_lo      <- 0.25         # Sensitive = bottom 25% AUC
auc_hi      <- 0.75         # Resistant = top 25% AUC
alpha_fdr   <- 0.05
logFC_cut   <- 0.58         # ≈1.5-fold (|log2FC|)
min_group   <- 5
label_top_n <- 10           # fallback labels if none of TOP_GENES are significant
png_res     <- 300

# EXACT drug names as in your GDSC file
focus_drugs <- c(
  "AZD8186","Avagacestat","BMS-536924","Afatinib","Foretinib",
  "AZD4547","AZD7762","AZD6738","Gemcitabine",
  "Docetaxel","Camptothecin","5-Fluorouracil"
)

# -------------------------
# HELPERS
# -------------------------
normalize_name <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]", "", x)
  toupper(x)
}
first_non_na <- function(...) {
  vals <- c(...)
  vals[which(!is.na(vals))[1]]
}

# -------------------------
# TOP GENES (for labeling)
# -------------------------
text_lines <- readLines(top_genes_file, warn = FALSE)
extract_symbol <- function(line) {
  toks <- unlist(strsplit(trimws(line), "\\s+"))
  toks[length(toks)]
}
TOP_GENES <- unique(na.omit(sapply(text_lines, extract_symbol)))

# -------------------------
# MODEL MAP (for harmonized names)
# -------------------------
model_df <- read_csv(model_path, show_col_types = FALSE) %>%
  mutate(
    SN = if ("StrippedCellLineName" %in% names(.)) StrippedCellLineName else NA,
    CN = if ("CellLineName" %in% names(.))         CellLineName         else NA,
    NormalizedName = coalesce(normalize_name(SN), normalize_name(CN))
  ) %>%
  select(ModelID, NormalizedName)

# -------------------------
# GDSC (AUC matrix)
# -------------------------
gdsc <- read_csv(gdsc_path, show_col_types = FALSE) %>%
  mutate(
    A = if ("cellname_modelID" %in% names(.)) normalize_name(cellname_modelID) else NA,
    B = if ("final_cellname"   %in% names(.)) normalize_name(final_cellname)   else NA,
    C = if ("CellLineName"     %in% names(.)) normalize_name(CellLineName)     else NA
  ) %>%
  left_join(model_df, by = "ModelID") %>%
  mutate(HarmonizedName = pmap_chr(list(A, B, C, NormalizedName), first_non_na))

# -------------------------
# EXPRESSION
# -------------------------
exp <- read_csv(exp_path, show_col_types = FALSE) %>%
  mutate(CellLineName_clean = normalize_name(CellLineName)) %>%
  distinct(CellLineName_clean, .keep_all = TRUE)

# -------------------------
# JOIN (expression + selected drugs)
# -------------------------
joined <- exp %>%
  inner_join(
    gdsc %>% select(HarmonizedName, all_of(focus_drugs)),
    by = c("CellLineName_clean" = "HarmonizedName")
  )

expr_cols <- setdiff(names(joined), c("CellLineName", "CellLineName_clean", focus_drugs))

# -------------------------
# LOOP: build per-drug volcanoes (without per-plot titles)
# -------------------------
plot_list <- list()
all_export <- list()   # to combine into one CSV later

for (drug in focus_drugs) {
  
  message("Processing: ", drug)
  
  auc <- joined[[drug]]
  qlo <- quantile(auc, auc_lo, na.rm = TRUE)
  qhi <- quantile(auc, auc_hi, na.rm = TRUE)
  
  joined2 <- joined %>%
    mutate(Group = ifelse(auc <= qlo, "Sensitive",
                          ifelse(auc >= qhi, "Resistant", NA))) %>%
    filter(!is.na(Group))
  
  # Skip if low counts
  if (length(unique(joined2$Group)) < 2 || min(table(joined2$Group)) < min_group) {
    message("  Skipped (insufficient samples)")
    next
  }
  
  # Reference level for direction
  joined2$Group <- factor(joined2$Group, levels = c("Sensitive", "Resistant"))
  
  # Expression matrix (genes x samples)
  E <- t(as.matrix(joined2[, expr_cols]))
  E <- E[apply(E, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
  if (nrow(E) == 0) {
    message("  Skipped (zero variance genes only)")
    next
  }
  
  # limma-trend on log-like expression
  design <- model.matrix(~ Group, data = joined2)
  fit <- lmFit(E, design)
  fit <- eBayes(fit, trend = TRUE)
  
  coef_name <- "GroupResistant"
  if (!(coef_name %in% colnames(fit$coefficients))) {
    message("  Skipped (missing contrast)")
    next
  }
  
  tt <- topTable(fit, coef = coef_name, number = Inf, p.value = 1, adjust.method = "BH")
  if (nrow(tt) == 0) {
    message("  Skipped (empty DE table)")
    next
  }
  
  tt$gene <- rownames(tt)
  tt <- tt %>%
    rename(log2FC = logFC, PValue = P.Value, FDR = adj.P.Val) %>%
    relocate(gene)
  
  tt$neglog10FDR <- -log10(pmax(tt$FDR, 1e-300))
  tt$sig <- (tt$FDR < alpha_fdr & abs(tt$log2FC) >= logFC_cut)
  
  # r = Pearson correlation between log2FC and -log10(FDR) for this panel
  r_value <- suppressWarnings(cor(tt$log2FC, tt$neglog10FDR, method = "pearson"))
  
  # Labels: TOP GENES first; else top-|log2FC| significant
  label_df <- tt %>% filter(sig & gene %in% TOP_GENES)
  if (nrow(label_df) == 0) {
    label_df <- tt %>%
      filter(sig) %>%
      slice_max(order_by = abs(log2FC), n = label_top_n, with_ties = FALSE)
  }
  
  # Export this drug's table
  out_csv <- paste0("VolcanoData_", drug, ".csv")
  write_csv(tt, out_csv)
  
  ###############################################################################
  # CREATE 3×4 GRID (each plot labeled with drug name)
  ###############################################################################
  
  if (length(plot_list) == 0) {
    stop("No plots generated. Check data.")
  }
  
  # Create grid WITH DRUG LABELS
  grid <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 4,
    nrow = 3,
    align = "hv",
    labels = names(plot_list),     # <-- DRUG NAMES APPEAR ON EACH PLOT
    label_size = 12,
    label_fontface = "bold",
    label_x = 0.02,
    label_y = 0.98,
    hjust = 0,
    vjust = 1
  )
  
  # ONE overall title
  overall_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Sensitive vs Resistant (bottom 25% vs top 25% AUC)  •  FDR < 0.05  •  Blue = Significant, Green = Non‑Significant",
      fontface = "bold",
      size = 16,
      x = 0.5, y = 0.98,
      hjust = 0.5, vjust = 1
    )
  
  final_page <- cowplot::plot_grid(
    overall_title,
    grid,
    ncol = 1,
    rel_heights = c(0.10, 1)
  ) 
  
  grid <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 4,
    nrow = 3,
    align = "hv",
    labels = names(plot_list),     # <-- THIS ADDS DRUG NAMES
    label_size = 14,
    label_fontface = "bold",
    label_x = 0.03,
    label_y = 0.97,
    hjust = 0,
    vjust = 1
  )
  
  
  # --- Save the 12-plot grid as PNG + PDF ---
  ggsave(
    filename = "VOLCANO_GRID_3x4.png",
    plot     = final_page,
    width    = 20,   # adjust if you want a smaller/larger figure
    height   = 15,
    dpi      = 300
  )
  
  ggsave(
    filename = "VOLCANO_GRID_3x4.pdf",
    plot     = final_page,
    width    = 20,
    height   = 15
  )
  
  ###############################################################################
  # HEATMAP OF DIFFERENTIAL EXPRESSION FOR TOP GENES USED IN VOLCANO PLOT
  ###############################################################################
  suppressPackageStartupMessages({
    library(pheatmap)
    library(dplyr)
    library(tidyr)
    library(readr)
  })
  
  # -------------------------
  # 1. Collect expression values ONLY for TOP GENES
  # -------------------------
  exp_mat <- exp %>%
    filter(CellLineName_clean %in% joined$CellLineName_clean) %>%
    select(CellLineName_clean, all_of(expr_cols))
  
  # Keep only TOP GENES that exist in expression matrix
  genes_for_heatmap <- intersect(TOP_GENES, expr_cols)
  
  if (length(genes_for_heatmap) == 0) {
    stop("ERROR: None of the TOP GENES are in the expression matrix.")
  }
  
  expr_sub <- joined %>%
    select(CellLineName_clean, all_of(genes_for_heatmap))
  
  expr_heat <- as.matrix(t(expr_sub[, genes_for_heatmap]))
  colnames(expr_heat) <- expr_sub$CellLineName_clean
  
  # Z‑score scale genes (rows)
  expr_heat_scaled <- t(scale(t(expr_heat)))
  
  # -------------------------
  # 2. Sensitivity/Resistance annotation (based on ANY drug)
  #    Here we use Afatinib just to provide colors for the heatmap annotation.
  # -------------------------
  drug_for_annotation <- "Afatinib"
  
  auc <- joined[[drug_for_annotation]]
  qlo <- quantile(auc, 0.25, na.rm = TRUE)
  qhi <- quantile(auc, 0.75, na.rm = TRUE)
  
  ann_df <- joined %>%
    mutate(Group = ifelse(auc <= qlo, "Sensitive",
                          ifelse(auc >= qhi, "Resistant", "Mid"))) %>%
    select(CellLineName_clean, Group)
  
  ann_df <- ann_df[match(colnames(expr_heat_scaled), ann_df$CellLineName_clean), , drop = FALSE]
  
  ann_colors <- list(
    Group = c(Sensitive = "#1f77b4", Resistant = "#d62728", Mid = "grey80")
  )
  
  # Force numeric matrix
  expr_heat_scaled <- apply(expr_heat_scaled, 2, as.numeric)
  expr_heat_scaled <- as.matrix(expr_heat_scaled)
  
  # Remove rows with all NA or NaN
  expr_heat_scaled <- expr_heat_scaled[rowSums(is.na(expr_heat_scaled)) < ncol(expr_heat_scaled), , drop = FALSE]
  
  # Remove columns with all NA (should not happen, but safe)
  expr_heat_scaled <- expr_heat_scaled[, colSums(is.na(expr_heat_scaled)) < nrow(expr_heat_scaled), drop = FALSE]
  
  test_numeric <- apply(expr_heat_scaled, 2, function(x) all(is.numeric(x)))
  test_numeric
  apply(expr_heat_scaled, 1, function(x) all(is.numeric(x)))
  bad_positions <- which(!apply(expr_heat_scaled, c(1,2), is.numeric), arr.ind = TRUE)
  bad_positions
  expr_heat_scaled[bad_positions]
  storage.mode(expr_heat_scaled)
  all(sapply(expr_heat_scaled, is.numeric))
  
  bad <- which(sapply(expr_heat_scaled, is.numeric) == FALSE, arr.ind = TRUE)
  bad
  expr_heat_scaled[bad]
  
  str(ann_df)
  table(ann_df$Group, useNA = "always")
  unique(ann_df$Group)
  
  # 1. Remove duplicate cell lines
  ann_unique <- ann_df %>%
    group_by(CellLineName_clean) %>%
    summarise(Group = first(Group)) %>%
    ungroup() %>%
    as.data.frame()
  
  # 2. Set rownames safely
  rownames(ann_unique) <- ann_unique$CellLineName_clean
  
  # 3. Align annotation to heatmap column order
  ann_aligned <- ann_unique[colnames(expr_heat_scaled), , drop = FALSE]
  
  # 4. Handle missing group labels (you had exactly ONE NA)
  ann_aligned$Group <- as.character(ann_aligned$Group)
  ann_aligned$Group[is.na(ann_aligned$Group)] <- "Unknown"
  
  # 5. Convert to the exact data.frame format pheatmap expects
  ann_aligned <- ann_aligned[, "Group", drop = FALSE]
  
  # Print first few names to inspect
  head(colnames(expr_heat_scaled))
  head(rownames(ann_aligned))
  
  # Which names are missing?
  setdiff(colnames(expr_heat_scaled), rownames(ann_aligned))
  
  # Which annotation names are extra?
  setdiff(rownames(ann_aligned), colnames(expr_heat_scaled)) 
  
  # Drop unmatched annotations like "KMH2.1"
  ann_aligned <- ann_aligned[colnames(expr_heat_scaled), , drop = FALSE]
  
  ###############################################################################
  # MAIN HEATMAP: TOP GENES (rows) × Samples (cols)
  # - Uses Z-scored expression per gene
  # - Column annotation = Group (Sensitive/Resistant/Mid/Unknown)
  # - Saves PNG + PDF
  ###############################################################################
  suppressPackageStartupMessages({
    library(pheatmap)
  })
  
  identical(rownames(ann_aligned), colnames(expr_heat_scaled))
 
  # Reuse your earlier normalizer (same as in your script)
  normalize_name <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("[^A-Za-z0-9]", "", x)  # remove non-alphanumerics
    toupper(x)
  }
  
  # Normalize both sides
  colnames(expr_heat_scaled) <- normalize_name(colnames(expr_heat_scaled))
  rownames(ann_aligned)      <- normalize_name(rownames(ann_aligned)) 
  
  # Compute the exact intersection
  common_samples <- intersect(colnames(expr_heat_scaled), rownames(ann_aligned))
  
  # Safety check: ensure we have samples left
  if (length(common_samples) == 0) {
    stop("No common samples between heatmap matrix and annotation after normalization.")
  }
  
  # Subset both objects to common_samples ONLY
  expr_heat_scaled <- expr_heat_scaled[, common_samples, drop = FALSE]
  ann_aligned      <- ann_aligned[common_samples, , drop = FALSE]
  
  # Order both by the same vector
  expr_heat_scaled <- expr_heat_scaled[, common_samples, drop = FALSE]
  ann_aligned      <- ann_aligned[common_samples, , drop = FALSE]
  
  # Make absolutely sure annotation is a plain data.frame, one column named "Group"
  ann_aligned <- as.data.frame(ann_aligned)
  if (!"Group" %in% colnames(ann_aligned)) {
    # If ann_aligned still contains CellLineName_clean, select Group explicitly
    ann_aligned <- ann_aligned[, "Group", drop = FALSE]
  } else {
    ann_aligned <- ann_aligned[, "Group", drop = FALSE]
  }
  
  # Final checks
  stopifnot(identical(rownames(ann_aligned), colnames(expr_heat_scaled)))
  stopifnot(is.data.frame(ann_aligned))
  stopifnot(storage.mode(expr_heat_scaled) %in% c("double", "numeric")) 
  
  # Ensure annotation colors contain all levels present
  if (!exists("ann_colors") || !"Group" %in% names(ann_colors)) {
    ann_colors <- list(Group = c(Sensitive="#1f77b4", Resistant="#d62728", Mid="grey80", Unknown="black"))
  } else {
    # Add Unknown if missing
    if (!"Unknown" %in% names(ann_colors$Group)) {
      ann_colors$Group <- c(ann_colors$Group, Unknown="black")
    }
  }
  
  pheatmap::pheatmap(
    expr_heat_scaled,
    annotation_col = ann_aligned,
    annotation_colors = ann_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    fontsize_col = 6,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    main = "TOP GENES Expression Heatmap (Z‑scored)"
  )
  
  ###############################################################################
  # FINAL HEATMAP EXPORT (PNG + PDF)
  ###############################################################################
  
  suppressPackageStartupMessages({
    library(pheatmap)
  })
  
  # --- Final alignment check ---
  common_samples <- intersect(colnames(expr_heat_scaled), rownames(ann_aligned))
  
  expr_heat_scaled <- expr_heat_scaled[, common_samples, drop = FALSE]
  ann_aligned      <- ann_aligned[common_samples, , drop = FALSE]
  
  stopifnot(identical(colnames(expr_heat_scaled), rownames(ann_aligned)))
  
  # --- Plot settings ---
  heat_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  
  ###############################################################################
  # TOP 50 GENES HEATMAP (labeled)
  # - Pick 50 genes by chosen method (variance across samples OR limma DE for a drug)
  # - Z-score per gene, label all rows, annotate columns by Group
  # - Saves PNG + PDF
  ###############################################################################
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(pheatmap)
  })
  
  # -------------------------
  # SETTINGS
  # -------------------------
  method <- "variance"        # options: "variance", "sig_drug", "abs_logFC"
  drug_for_top50 <- "Afatinib"  # used if method is sig_drug or abs_logFC
  
  # -------------------------
  # INPUTS assumed already available from your pipeline:
  #   exp, ann_df, TOP_GENES, joined, and (optionally) tt per drug if using sig_drug/abs_logFC
  # If you have per-drug tt tables saved as CSVs (e.g., VolcanoData_<drug>.csv), we can read them.
  # -------------------------
  
  # 1) Expression gene list from 'exp'
  expr_cols <- setdiff(names(exp), c("CellLineName","CellLineName_clean"))
  
  # 2) Candidate genes limited to TOP_GENES present in expression
  top_genes_in_matrix <- intersect(TOP_GENES, expr_cols)
  if (length(top_genes_in_matrix) < 50) {
    warning(sprintf("Only %d TOP_GENES found in expression. Will use all of them.", length(top_genes_in_matrix)))
  }
  
  # 3) Build expression matrix for TOP_GENES (genes x samples)
  expr_top <- exp %>%
    select(CellLineName_clean, all_of(top_genes_in_matrix)) %>%
    distinct(CellLineName_clean, .keep_all = TRUE) %>%
    column_to_rownames("CellLineName_clean") %>%
    t() %>%
    as.matrix()
  
  # Force numeric and drop zero-variance genes
  expr_top <- apply(expr_top, 2, function(x) as.numeric(as.character(x)))
  expr_top <- as.matrix(expr_top)
  expr_top <- expr_top[apply(expr_top, 1, function(v) var(v, na.rm=TRUE) > 0), , drop = FALSE]
  
  # -------------------------
  # 4) Choose the TOP 50 genes
  # -------------------------
  pick_top50 <- function(expr_mat, method = "variance", drug = NULL) {
    if (method == "variance") {
      # Top 50 most variable genes across samples
      v <- apply(expr_mat, 1, function(v) stats::var(v, na.rm = TRUE))
      genes50 <- names(sort(v, decreasing = TRUE))[seq_len(min(50, length(v)))]
      return(genes50)
    }
    # Methods based on a chosen drug use limma table 'tt'
    if (is.null(drug)) stop("Please set 'drug_for_top50' for methods 'sig_drug' or 'abs_logFC'.")
    tt_path <- paste0("VolcanoData_", drug, ".csv")
    if (!file.exists(tt_path)) {
      stop("Could not find ", tt_path, 
           ". Re-run the volcano script that writes VolcanoData_<drug>.csv, or switch method='variance'.")
    }
    tt <- suppressMessages(readr::read_csv(tt_path, show_col_types = FALSE))
    # Ensure gene column naming
    if (!"gene" %in% names(tt)) {
      if ("X1" %in% names(tt)) names(tt)[names(tt)=="X1"] <- "gene"
      if (!"gene" %in% names(tt)) stop("Limma table must have a 'gene' column.")
    }
    # Keep only genes present in expression
    tt <- tt %>% filter(gene %in% rownames(expr_mat))
    if (nrow(tt) == 0) stop("No overlap between limma table and expression genes.")
    if (method == "sig_drug") {
      # Significant by FDR & |log2FC|, then strongest by -log10(FDR)
      if (!all(c("FDR","log2FC") %in% names(tt))) stop("Limma table missing FDR/log2FC.")
      tt <- tt %>%
        mutate(sig = (FDR < 0.05 & abs(log2FC) >= 0.58),
               score = -log10(pmax(FDR, 1e-300)))
      sig_tt <- tt %>% filter(sig)
      if (nrow(sig_tt) == 0) {
        warning("No significant genes at FDR<0.05 & |log2FC|>=0.58; falling back to top by score.")
        sig_tt <- tt %>% mutate(score = -log10(pmax(FDR, 1e-300)))
      }
      genes50 <- sig_tt %>%
        arrange(desc(score)) %>%
        slice_head(n = 50) %>%
        pull(gene)
      return(genes50)
    } else if (method == "abs_logFC") {
      if (!"log2FC" %in% names(tt)) stop("Limma table missing log2FC.")
      genes50 <- tt %>%
        arrange(desc(abs(log2FC))) %>%
        slice_head(n = 50) %>%
        pull(gene)
      return(genes50)
    } else {
      stop("Unknown method: ", method)
    }
  }
  
  genes50 <- pick_top50(expr_top, method = method, drug = drug_for_top50)
  
  # Subset to the selected 50 (or fewer if not enough)
  expr_50 <- expr_top[intersect(genes50, rownames(expr_top)), , drop = FALSE]
  
  # -------------------------
  # 5) Z-score rows and drop all-NA rows after scaling
  # -------------------------
  expr_50_z <- t(scale(t(expr_50)))
  expr_50_z <- expr_50_z[rowSums(is.na(expr_50_z)) < ncol(expr_50_z), , drop = FALSE]
  
  # -------------------------
  # 6) Build & align column annotation (Group)
  # -------------------------
  ann_unique <- ann_df %>%
    group_by(CellLineName_clean) %>%
    summarise(Group = first(Group)) %>%
    ungroup() %>%
    as.data.frame()
  rownames(ann_unique) <- ann_unique$CellLineName_clean
  
  common_samples <- intersect(colnames(expr_50_z), rownames(ann_unique))
  expr_50_z  <- expr_50_z[, common_samples, drop = FALSE]
  ann_aligned <- ann_unique[common_samples, , drop = FALSE]
  
  ann_aligned$Group <- as.character(ann_aligned$Group)
  ann_aligned$Group[is.na(ann_aligned$Group)] <- "Unknown"
  ann_aligned <- ann_aligned[, "Group", drop = FALSE]
  
  
  # ============================================================
  #  TOP 50 GENES × 12 DRUGS — SPEARMAN CORRELATION HEATMAP (v2)
  #  Files required in working directory:
  #    - Expression_clean_harmonized.csv
  #    - GDSC_AUC_with_cellname_modelID.csv
  #    - TOP GENES.txt
  # ============================================================
  
  # ---- 0) Install / load packages ----
  ensure_pkg <- function(pkg, bioc = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg)
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  ensure_pkg("tidyverse")
  ensure_pkg("stringr")
  ensure_pkg("circlize")                 # for ComplexHeatmap colors
  ok_ch <- tryCatch({ ensure_pkg("ComplexHeatmap", bioc = TRUE); TRUE }, error = function(e) FALSE)
  if (!ok_ch) ensure_pkg("pheatmap")     # fallback
  
  # ---- 1) Your 12 drugs as typed (typos will be handled) ----
  target_input <- c(
    "AZD8186", "Avagacestat", "BMS356924", "Afatanib",
    "foretinib", "AZD4547", "AZD7762", "AZD6738",
    "Gemecitabine", "docetaxel", "campothecin", "5 flouracil"
  )
  
  # ---- 2) Load data ----
  expr_raw <- readr::read_csv("Expression_clean_harmonized.csv", show_col_types = FALSE)
  drug_raw <- readr::read_csv("GDSC_AUC_with_cellname_modelID.csv", show_col_types = FALSE)
  top_raw  <- readr::read_lines("TOP GENES.txt")
  
  # ---- 3) Parse top-genes (strip ranks; keep symbols only) ----
  top_genes <- top_raw %>%
    stringr::str_replace("^\\s*\\d+\\s+", "") %>%
    stringr::str_trim() %>%
    { .[. != ""] }
  
  # ---- 4) Prepare expression: collapse duplicates; transpose to genes × samples ----
  stopifnot(all(c("CellLineName", "CellLineName_clean") %in% colnames(expr_raw)))  # file layout check  [1](https://mylifeunisaac-my.sharepoint.com/personal/22994327_mylife_unisa_ac_za/_layouts/15/Doc.aspx?sourcedoc=%7B8D23F68B-179F-4596-BF77-6D77FE8BA3B0%7D&file=Expression_clean_harmonized.csv&action=default&mobileredirect=true)
  
  gene_cols <- setdiff(colnames(expr_raw), c("CellLineName", "CellLineName_clean"))
  expr_tmp <- expr_raw %>%
    mutate(sample_id = dplyr::if_else(is.na(CellLineName_clean) | CellLineName_clean == "",
                                      CellLineName, CellLineName_clean)) %>%
    select(sample_id, dplyr::all_of(gene_cols))
  
  dups <- expr_tmp %>% count(sample_id) %>% filter(n > 1)
  if (nrow(dups) > 0) {
    message("Aggregating duplicated sample IDs by mean:\n",
            paste0("  - ", dups$sample_id, " (", dups$n, " rows)", collapse = "\n"))
  }
  expr_collapsed <- expr_tmp %>%
    group_by(sample_id) %>%
    summarise(across(dplyr::all_of(gene_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  expr_samples_mat <- expr_collapsed %>% column_to_rownames("sample_id") %>% as.data.frame()
  expr <- t(as.matrix(expr_samples_mat))  # genes × samples
  
  # ---- 5) Robust drug-name mapping (normalize + synonyms + approximate) ----
  normalize_name <- function(x) tolower(gsub("[^a-z0-9]+", "", x))
  synonyms <- c(  # fix common variants to canonical GDSC spellings  [2](https://mylifeunisaac-my.sharepoint.com/personal/22994327_mylife_unisa_ac_za/_layouts/15/Doc.aspx?sourcedoc=%7B57CAC09E-7602-4DD8-B89E-89879DEAD321%7D&file=GDSC_AUC_with_cellname_modelID.csv&action=default&mobileredirect=true)
    "bms356924"   = "bms536924",
    "afatanib"    = "afatinib",
    "gemecitabine"= "gemcitabine",
    "campothecin" = "camptothecin",
    "5flouracil"  = "5fluorouracil"
  )
  data_cols <- colnames(drug_raw)
  data_norm <- normalize_name(data_cols)
  map_one <- function(req) {
    req_norm <- normalize_name(req)
    if (req_norm %in% names(synonyms)) req_norm <- synonyms[[req_norm]]
    exact_ix   <- which(data_norm == req_norm);            if (length(exact_ix) > 0) return(data_cols[exact_ix[1]])
    starts_ix  <- which(startsWith(data_norm, req_norm));  if (length(starts_ix) > 0) return(data_cols[starts_ix[1]])
    contains_ix<- which(grepl(req_norm, data_norm, fixed = TRUE)); if (length(contains_ix) > 0) return(data_cols[contains_ix[1]])
    approx_ix  <- agrep(req_norm, data_norm, max.distance = 0.2);  if (length(approx_ix) > 0) return(data_cols[approx_ix[1]])
    return(NA_character_)
  }
  mapping <- tibble::tibble(
    requested = target_input,
    matched_column = vapply(target_input, map_one, character(1))
  )
  cat("\n== Requested → Matched GDSC columns ==\n"); print(mapping, n = nrow(mapping))
  mapping_ok <- mapping %>% filter(!is.na(matched_column))
  
  # ---- 6) Build drug_use (safe select + rename) ----
  stopifnot("final_cellname" %in% colnames(drug_raw))  # column exists in your file  [2](https://mylifeunisaac-my.sharepoint.com/personal/22994327_mylife_unisa_ac_za/_layouts/15/Doc.aspx?sourcedoc=%7B57CAC09E-7602-4DD8-B89E-89879DEAD321%7D&file=GDSC_AUC_with_cellname_modelID.csv&action=default&mobileredirect=true)
  drug_sel <- dplyr::select(drug_raw, final_cellname, dplyr::any_of(mapping_ok$matched_column))
  old_names <- intersect(mapping_ok$matched_column, colnames(drug_sel))
  new_names <- mapping_ok$requested[match(old_names, mapping_ok$matched_column)]
  name_map  <- setNames(new_names, old_names)
  for (old in names(name_map)) names(drug_sel)[names(drug_sel) == old] <- name_map[[old]]
  
  # ---- 7) ***AUTO-SELECT BEST SAMPLE-ID COLUMN*** and align ----
  # Your AUC file provides several ID columns; we pick the one that maximizes overlap with expression.  [2](https://mylifeunisaac-my.sharepoint.com/personal/22994327_mylife_unisa_ac_za/_layouts/15/Doc.aspx?sourcedoc=%7B57CAC09E-7602-4DD8-B89E-89879DEAD321%7D&file=GDSC_AUC_with_cellname_modelID.csv&action=default&mobileredirect=true)
  id_cols_in_drug <- c("final_cellname", "CellLineName", "GDSC_CellLineName", "cellname_modelID")
  id_cols_in_drug <- id_cols_in_drug[id_cols_in_drug %in% colnames(drug_raw)]  # keep those that exist
  
  # normalized expression IDs (from expr column names)
  expr_ids_norm <- normalize_name(colnames(expr))
  
  # compute overlaps for each candidate drug ID column
  overlaps <- sapply(id_cols_in_drug, function(col) {
    ids <- normalize_name(drug_raw[[col]])
    length(intersect(expr_ids_norm, ids))
  })
  best_col <- id_cols_in_drug[which.max(overlaps)]
  cat("\nBest matching ID column in AUC file:", best_col, "with overlap =", max(overlaps), "\n")
  
  # attach a normalized match_id from the chosen column
  drug_sel$match_id <- normalize_name(drug_raw[[best_col]])
  drug_sel <- dplyr::filter(drug_sel, !is.na(match_id) & match_id != "")
  
  # ---------------------------
  # 6) Build drug_use (safe select + rename)  ——  REPLACEMENT
  # ---------------------------
  stopifnot("final_cellname" %in% colnames(drug_raw))  # exists in your file  [1](https://mylifeunisaac-my.sharepoint.com/personal/22994327_mylife_unisa_ac_za/_layouts/15/Doc.aspx?sourcedoc=%7B57CAC09E-7602-4DD8-B89E-89879DEAD321%7D&file=GDSC_AUC_with_cellname_modelID.csv&action=default&mobileredirect=true)
  
  # Select the matched drug columns + a few candidate ID columns for alignment
  candidate_ids <- c("final_cellname", "CellLineName", "cellname_modelID", "GDSC_CellLineName")
  have_ids      <- candidate_ids[candidate_ids %in% colnames(drug_raw)]
  drug_sel <- dplyr::select(drug_raw, dplyr::any_of(have_ids), dplyr::any_of(mapping_ok$matched_column))
  
  # Rename matched drug columns to your requested names (outside the pipe)
  old_names <- intersect(mapping_ok$matched_column, colnames(drug_sel))
  new_names <- mapping_ok$requested[match(old_names, mapping_ok$matched_column)]
  name_map  <- setNames(new_names, old_names)
  for (old in names(name_map)) {
    names(drug_sel)[names(drug_sel) == old] <- name_map[[old]]
  }
  
  cat("\nDrug columns selected (post-rename):\n")
  print(setdiff(colnames(drug_sel), have_ids))
  
  # ---------------------------
  # 7) Choose ID column and align ——  REPLACEMENT
  # ---------------------------
  
  # Conservative cleaner: keep letters/digits, lower-case, collapse spaces/hyphens to one underscore
  clean_id <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub("[[:space:]]+", "_", x)
    x <- gsub("[-]+", "_", x)
    x <- gsub("[^a-z0-9_]+", "", x)  # keep letters, digits, underscore
    x
  }
  
  # Build *expression* ID vector to compare against
  expr_ids_clean <- clean_id(colnames(expr))
  
  # Try candidate ID columns in a sensible order and pick the first with decent overlap (>= 5).
  best_col <- NA_character_
  best_overlap <- -1L
  for (col in have_ids) {
    # Some ID columns might be numeric; coerce to character safely
    ids <- as.character(drug_sel[[col]])
    ids_clean <- clean_id(ids)
    ov <- length(intersect(expr_ids_clean, ids_clean))
    cat(sprintf("Candidate AUC ID column '%s' → overlap = %d\n", col, ov))
    if (ov > best_overlap) {
      best_overlap <- ov
      best_col <- col
    }
    if (ov >= 5) break
  }
  cat(sprintf("\nChosen AUC ID column: %s (overlap=%d)\n", best_col, best_overlap))
  
  # Keep rows with a non-empty chosen ID
  ids_clean <- clean_id(as.character(drug_sel[[best_col]]))
  drug_sel2 <- drug_sel %>%
    mutate(.match_id = ids_clean) %>%
    filter(!is.na(.match_id) & .match_id != "")
  
  # Aggregate duplicates by .match_id (mean AUC across replicates)
  drug_agg <- drug_sel2 %>%
    group_by(.match_id) %>%
    summarise(across(.cols = where(is.numeric), .fns = ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  # Convert to data.frame and set rownames safely
  drug_use <- as.data.frame(drug_agg)
  rownames(drug_use) <- drug_use$.match_id
  drug_use$.match_id <- NULL
  
  # Final list of drug columns (should be your 12, or subset if some didn’t map)
  cat("\nFinal drug columns in drug_use:\n")
  print(colnames(drug_use))
  
  # ---------------------------
  # 8) Align samples  ——  REPLACEMENT
  # ---------------------------
  # Clean expression sample IDs the same way (already in expr_ids_clean)
  colnames(expr) <- expr_ids_clean
  
  common_samples <- intersect(colnames(expr), rownames(drug_use))
  cat("\n# Common samples between expression and AUC: ", length(common_samples), "\n")
  
  if (length(common_samples) < 5) {
    message(
      "\nOverlap is still small. Quick hints:\n",
      " - Inspect a few expression IDs vs AUC IDs:\n",
      "   head(colnames(expr)); head(rownames(drug_use))\n",
      " - If needed, force a specific ID column by setting best_col manually to one of:\n",
      sprintf("   %s\n", paste(have_ids, collapse = ", "))
    )
    # You can also print a few nearest suggestions:
    ex_ids <- setdiff(colnames(expr), rownames(drug_use))[1:min(10, length(setdiff(colnames(expr), rownames(drug_use))))]
    nearest <- function(x, pool) pool[which.min(adist(x, pool))]
    if (length(ex_ids) > 0) {
      cat("\nNearest AUC IDs for a few expression IDs:\n")
      for (x in ex_ids) cat("  ", x, "  →  ", nearest(x, rownames(drug_use)), "\n")
    }
    stop("Too few overlapping samples after ID harmonization. See hints above.")
  }
  
  expr_aligned <- expr[, common_samples, drop = FALSE]
  drug_aligned <- as.matrix(drug_use[common_samples, , drop = FALSE])
  
  # normalize expression column IDs to the same scheme
  colnames(expr) <- expr_ids_norm
  
  # compute overlap and show a few nearest suggestions if overlap is small
  common_samples <- intersect(colnames(expr), rownames(drug_use))
  cat("\n# Common samples between expression and AUC: ", length(common_samples), "\n")
  
  if (length(common_samples) < 5) {
    message("\nOverlap is small. Showing a few nearest suggestions for 10 expression IDs:")
    ex_ids <- setdiff(colnames(expr), rownames(drug_use))
    ex_ids <- ex_ids[1:min(10, length(ex_ids))]
    nearest <- function(x, pool) pool[which.min(adist(x, pool))]
    for (x in ex_ids) cat("  ", x, "  →  ", nearest(x, rownames(drug_use)), "\n")
    stop("Too few overlapping samples after ID harmonization. Inspect suggestions and naming.")
  }
  
  # Now subset aligned matrices
  expr_aligned <- expr[, common_samples, drop = FALSE]
  drug_aligned <- as.matrix(drug_use[common_samples, , drop = FALSE])
  
  # ---- 8) Keep only top genes present in expression ----
  top_available <- intersect(top_genes, rownames(expr_aligned))
  if (length(top_available) < 10) warning("Fewer than 10 top genes found; check symbols in TOP GENES.txt.")
  expr_top <- expr_aligned[top_available, , drop = FALSE]
  
  # DGE heatmap (Sensitive vs Resistant) for a specified gene list
  # - drug_name: one of colnames(auc_aln)
  # - genes_to_plot: character vector of symbols (e.g., top_genes)
  # - q: quantile cutoff for S vs R
  run_dge_heatmap_for_genes <- function(drug_name, genes_to_plot, q = 0.25, min_per_group = 3, topN = NULL) {
    stopifnot(drug_name %in% colnames(auc_aln))
    suppressPackageStartupMessages({
      library(limma)
      if (ok_ch && requireNamespace("ComplexHeatmap", quietly = TRUE)) library(ComplexHeatmap)
    })
    
    # ---- define Sensitive / Resistant by AUC quantiles ----
    y <- as.numeric(auc_aln[, drug_name]); names(y) <- rownames(auc_aln)
    lo <- quantile(y, probs = q, na.rm = TRUE)
    hi <- quantile(y, probs = 1 - q, na.rm = TRUE)
    sens_ids <- names(y)[y <= lo]
    resi_ids <- names(y)[y >= hi]
    keep_ids <- c(sens_ids, resi_ids)
    if (length(sens_ids) < min_per_group || length(resi_ids) < min_per_group) {
      message(sprintf("[Skip %s] Too few samples: S=%d, R=%d", drug_name, length(sens_ids), length(resi_ids)))
      return(invisible(NULL))
    }
    
    # ---- subset expression to the kept samples ----
    X <- expr_aln[, keep_ids, drop = FALSE]     # genes × kept samples
    group <- factor(ifelse(colnames(X) %in% resi_ids, "Resistant", "Sensitive"),
                    levels = c("Sensitive", "Resistant"))
    
    # ---- limma DGE to obtain stats (optional but proper DGE) ----
    design <- model.matrix(~ group)
    fit <- lmFit(X, design)
    fit <- eBayes(fit, robust = TRUE)
    tt  <- topTable(fit, coef = 2, n = Inf, sort.by = "P")  # Resistant vs Sensitive
    
    # ---- choose which genes to show ----
    # 1) user-specified list, intersect with what’s present
    g_list <- intersect(genes_to_plot, rownames(X))
    # 2) optionally override with topN DE genes by FDR if requested
    if (!is.null(topN)) {
      top_dge <- rownames(tt)[order(tt$adj.P.Val, -abs(tt$logFC))][seq_len(min(topN, nrow(tt)))]
      g_list <- intersect(top_dge, rownames(X))
    }
    if (length(g_list) < 2) {
      message(sprintf("[Skip %s] <2 genes available after intersection; check symbols.", drug_name))
      return(invisible(NULL))
    }
    
    # ---- build Z-scored heatmap matrix (genes × kept samples) ----
    Z <- t(scale(t(X[g_list, , drop = FALSE])))
    Z[is.na(Z)] <- 0
    annot <- data.frame(Phenotype = group, row.names = colnames(Z))
    
    # ---- plot & save (PNG + PDF) ----
    base <- gsub("[^A-Za-z0-9]+", "_", drug_name)
    
    if (ok_ch && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      ha <- HeatmapAnnotation(
        df  = annot,
        col = list(Phenotype = c(Sensitive = "#1f78b4", Resistant = "#e31a1c"))
      )
      
      ht <- Heatmap(
        Z, name = "Z-score",
        top_annotation = ha,
        col = circlize::colorRamp2(c(-2, 0, 2), c("#4575b4", "white", "#d73027")),
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 8),
        column_title = sprintf("%s — DGE heatmap (Top %d)", drug_name, topN),
        row_title = "DE genes"
      )  # <- CLOSE Heatmap() here
      
      # Save PNG
      png(sprintf("DGE_Heatmap_%s.png", base), width = 2000, height = 2400, res = 300)
      ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
      dev.off()
      
      # Save PDF
      pdf(sprintf("DGE_Heatmap_%s.pdf", base), width = 8.5, height = 10.5)
      ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
      dev.off()
      
    } else {
      # pheatmap fallback
      pheatmap::pheatmap(
        Z, annotation_col = annot, show_colnames = FALSE,
        color = colorRampPalette(c("#4575b4","white","#d73027"))(100),
        fontsize_row = 8, clustering_method = "average",
        filename = sprintf("DGE_Heatmap_%s.png", base), width = 10, height = 12
      )
      grDevices::pdf(sprintf("DGE_Heatmap_%s.pdf", base), width = 8.5, height = 10.5)
      pheatmap::pheatmap(
        Z, annotation_col = annot, show_colnames = FALSE,
        color = colorRampPalette(c("#4575b4","white","#d73027"))(100),
        fontsize_row = 8, clustering_method = "average"
      )
      dev.off()
    }