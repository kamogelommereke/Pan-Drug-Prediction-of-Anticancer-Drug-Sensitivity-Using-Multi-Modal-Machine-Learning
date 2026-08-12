setwd(setwd("C:/DESKTOP/BIOINFORMATICS VALIDATION"))

# ===========================
# Validate DL-predicted genes:
# Expression–Drug correlation (heatmap + volcano)
# Focused on genes listed in 'TOP GENES.txt'
# ===========================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(scales)
  library(tibble)   # column_to_rownames
})

# ---------------------------
# User parameters (edit here)
# ---------------------------

# Input paths
exp_path_raw    <- "Expression_Transposed_CellNames_RENAMED.csv"
exp_path_clean  <- "Expression_clean_harmonized.csv"  # if available
gdsc_path       <- "GDSC_AUC_with_cellname_modelID.csv"
model_path      <- "Model.csv"
top_genes_path  <- "TOP GENES.txt"                    # your uploaded file

# Use already harmonized expression file if present?
use_harmonized_expression <- TRUE

# Correlation options
corr_method   <- "spearman"   # "pearson" or "spearman"
min_samples   <- 30           # min overlapping cell lines per gene–drug correlation
alpha_fdr     <- 0.05         # FDR threshold for volcano/heatmap highlights
effect_cut    <- 0.25         # |r| threshold for highlights

# Drug selection for heatmap
drugs_of_interest <- NULL     # e.g., c("Cisplatin","Paclitaxel","Doxorubicin"); NULL -> auto-pick by coverage
max_genes_heatmap <- 300      # cap for heatmap rows (your gene list is smaller, so usually not hit)

# Volcano plot settings
volcano_drug  <- NULL         # set a single drug; if NULL, uses first in 'drugs_of_interest'
label_top_n   <- 15           # label top-N genes by |r| among significant

# AUC handling: default AUC lower = more sensitive.
# If you prefer "higher values = more sensitive" (for intuition in plots), set invert_auc = TRUE.
invert_auc    <- FALSE

# File output
png_width  <- 2200; png_height <- 1800; png_res <- 250
pdf_width  <- 10;   pdf_height <- 8

# ---------------------------
# Helpers
# ---------------------------

normalize_name <- function(x) {
  x <- trimws(as.character(x))
  x <- ifelse(is.na(x) | x == "NA", NA_character_, x)
  ifelse(is.na(x), NA_character_, str_replace_all(x, "[^0-9A-Za-z]+", ""))
}

first_non_na <- function(...) {
  vals <- c(...)
  if (all(is.na(vals))) NA_character_ else vals[which(!is.na(vals))[1]]
}

# Compute correlation vector of all genes vs one drug efficiently
cor_gene_vs_drug <- function(expr_mat, drug_vec, method = "spearman", min_n = 30) {
  ok_drug <- !is.na(drug_vec)
  if (sum(ok_drug) < min_n) {
    return(tibble(gene = colnames(expr_mat), r = NA_real_, p = NA_real_, n = 0L))
  }
  X <- expr_mat[ok_drug, , drop = FALSE]
  y <- drug_vec[ok_drug]
  
  n_vec <- colSums(!is.na(X) & !is.na(y))
  keep  <- n_vec >= min_n
  res <- tibble(gene = colnames(expr_mat), r = NA_real_, p = NA_real_, n = n_vec)
  if (!any(keep)) return(res)
  
  Xk <- X[, keep, drop = FALSE]
  r_vec <- suppressWarnings(cor(Xk, y, use = "pairwise.complete.obs", method = method))
  nk   <- n_vec[keep]
  r2   <- pmin(pmax(r_vec, -0.999999), 0.999999)
  tval <- r2 * sqrt((nk - 2) / (1 - r2^2))
  pval <- 2 * pt(abs(tval), df = pmax(nk - 2, 1), lower.tail = FALSE)
  
  res$r[keep] <- as.numeric(r_vec)
  res$p[keep] <- as.numeric(pval)
  res
}

# ---------------------------
# Read TOP GENES (predictions)
# ---------------------------
stopifnot(file.exists(top_genes_path))
top_genes_raw <- readLines(top_genes_path, warn = FALSE)

# your file appears like "1 TMEM45B", "2 GMFG", ...
# extract last token on each line as gene symbol (robust to tab/space)
extract_symbol <- function(line) {
  toks <- strsplit(trimws(line), "\\s+")[[1]]
  toks[length(toks)]
}
library(readr)
genes_tbl <- read_table("TOP GENES.txt", col_names = c("idx", "symbol"),
                        col_types = cols(idx = col_character(), symbol = col_character()),
                        comment = "#")  # skip any comment lines
top_genes <- unique(na.omit(genes_tbl$symbol))
# ---------------------------
# Load and harmonize
# ---------------------------

# 1) ModelID -> normalized cell name
model_df <- read_csv(model_path, show_col_types = FALSE)
stopifnot("ModelID" %in% names(model_df))
has_stripped <- "StrippedCellLineName" %in% names(model_df)
has_cellname <- "CellLineName" %in% names(model_df)

model_map <- model_df %>%
  mutate(NormalizedName = dplyr::case_when(
    has_stripped ~ normalize_name(StrippedCellLineName),
    has_cellname ~ normalize_name(CellLineName),
    TRUE ~ NA_character_
  )) %>%
  select(ModelID, NormalizedName) %>%
  distinct()

# 2) GDSC: build HarmonizedName & detect drug columns
gdsc_full <- read_csv(gdsc_path, show_col_types = FALSE, guess_max = 1e5) %>%
  mutate(
    cellname_modelID_clean = if ("cellname_modelID" %in% names(.)) normalize_name(cellname_modelID) else NA_character_,
    final_cellname_clean   = if ("final_cellname"   %in% names(.)) normalize_name(final_cellname)   else NA_character_,
    CellLineName_clean     = if ("CellLineName"     %in% names(.)) normalize_name(CellLineName)     else NA_character_
  )

if ("ModelID" %in% names(gdsc_full)) {
  gdsc_full <- gdsc_full %>%
    left_join(model_map, by = "ModelID") %>%
    rename(from_model_map_clean = NormalizedName)
} else {
  gdsc_full <- gdsc_full %>% mutate(from_model_map_clean = NA_character_)
}

priority_cols <- c("cellname_modelID_clean", "final_cellname_clean", "CellLineName_clean", "from_model_map_clean")
present_cols  <- priority_cols[priority_cols %in% names(gdsc_full)]
stopifnot(length(present_cols) > 0)

gdsc_full <- gdsc_full %>%
  mutate(HarmonizedName = pmap_chr(across(all_of(present_cols)), first_non_na))

name_cols <- c("GDSC_CellLineName","ModelID","CellLineName","final_cellname","cellname_modelID",
               "cellname_modelID_clean","final_cellname_clean","CellLineName_clean","from_model_map_clean",
               "HarmonizedName")
drug_cols <- setdiff(names(gdsc_full)[sapply(gdsc_full, is.numeric)], name_cols)

if (invert_auc) {
  gdsc_full[drug_cols] <- lapply(gdsc_full[drug_cols], function(z) ifelse(is.na(z), NA_real_, -z))
}

# 3) Expression
if (use_harmonized_expression && file.exists(exp_path_clean)) {
  exp_df <- read_csv(exp_path_clean, show_col_types = FALSE)
  if (!"CellLineName_clean" %in% names(exp_df) && "CellLineName" %in% names(exp_df)) {
    exp_df <- exp_df %>% mutate(CellLineName_clean = normalize_name(CellLineName))
  }
} else {
  exp_df <- read_csv(exp_path_raw, show_col_types = FALSE)
  stopifnot("CellLineName" %in% names(exp_df))
  exp_df <- exp_df %>% mutate(CellLineName_clean = normalize_name(CellLineName))
}

# de-duplicate cell lines
exp_df <- exp_df %>%
  filter(!is.na(CellLineName_clean)) %>%
  distinct(CellLineName_clean, .keep_all = TRUE)

# Diagnostics: overlaps
exp_names  <- exp_df %>% distinct(CellLineName_clean) %>% pull()
gdsc_names <- gdsc_full %>% distinct(HarmonizedName) %>% pull()
message("Overlap cell lines: ", length(intersect(exp_names, gdsc_names)))

# ---------------------------
# Join and build matrices
# ---------------------------
joined <- exp_df %>%
  inner_join(
    gdsc_full %>% select(HarmonizedName, all_of(drug_cols)),
    by = c("CellLineName_clean" = "HarmonizedName")
  )

if (nrow(joined) < min_samples) {
  stop("Too few overlapping cell lines after join. Consider lowering min_samples or checking name mapping.")
}

expr_name_cols <- c("CellLineName", "CellLineName_clean")
expr_candidates <- setdiff(names(joined), c(expr_name_cols, drug_cols))

# keep only numeric gene columns
expr_candidates <- expr_candidates[sapply(joined[expr_candidates], is.numeric)]

# ---------------------------
# Restrict to your predicted genes
# ---------------------------
# Match gene symbols case-sensitively first, then case-insensitive fallback
found_genes <- intersect(expr_candidates, top_genes)

if (length(found_genes) == 0) {
  # case-insensitive map
  map_ci <- setNames(expr_candidates, toupper(expr_candidates))
  want   <- unique(toupper(top_genes))
  found_ci <- unique(na.omit(map_ci[intersect(names(map_ci), want)]))
  found_genes <- found_ci
}

if (length(found_genes) == 0) {
  stop("None of the TOP GENES were found as numeric columns in your expression matrix. Check gene symbols / column names.")
}

if (length(found_genes) > max_genes_heatmap) {
  message("Your gene list (", length(found_genes), ") exceeds max_genes_heatmap (",
          max_genes_heatmap, "). Capping to first ", max_genes_heatmap, ".")
  found_genes <- found_genes[seq_len(max_genes_heatmap)]
}

message("Found ", length(found_genes), " genes from your DL predictions in the expression matrix.")

expr_mat <- as.matrix(joined[found_genes])
drug_mat <- joined[drug_cols]

# ---------------------------
# Select drugs for heatmap
# ---------------------------
if (is.null(drugs_of_interest)) {
  comp_counts  <- sapply(drug_mat, function(z) sum(!is.na(z)))
  drugs_of_interest <- names(sort(comp_counts, decreasing = TRUE))[1:min(12, length(comp_counts))]
  if (length(drugs_of_interest) == 0) stop("No numeric drug columns found in GDSC matrix.")
  message("Auto-selected drugs by coverage: ", paste(drugs_of_interest, collapse = ", "))
} else {
  drugs_of_interest <- intersect(drugs_of_interest, drug_cols)
  if (length(drugs_of_interest) == 0) stop("Specified drugs not found among GDSC columns.")
}

# ---------------------------
# Correlations (only your genes x selected drugs)
# ---------------------------
all_cor_long <- vector("list", length(drugs_of_interest))
names(all_cor_long) <- drugs_of_interest

for (d in drugs_of_interest) {
  dv <- joined[[d]]
  cor_res <- cor_gene_vs_drug(expr_mat, dv, method = corr_method, min_n = min_samples) %>%
    mutate(drug = d)
  all_cor_long[[d]] <- cor_res
}

cor_long <- bind_rows(all_cor_long) %>%
  group_by(drug) %>%
  mutate(FDR = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(direction = case_when(
    is.na(r) ~ "NA",
    r >= 0   ~ "positive",
    TRUE     ~ "negative"
  ))

write_csv(cor_long, "correlation_long_TOP_GENES_only.csv")

# ---------------------------
# Heatmap (genes = your list)
# ---------------------------
heat_df <- cor_long %>%
  filter(drug %in% drugs_of_interest, !is.na(r)) %>%
  select(drug, gene, r) %>%
  distinct() %>%
  pivot_wider(names_from = drug, values_from = r)

write_csv(as.data.frame(heat_df), "correlation_wide_heatmap_matrix_TOP_GENES.csv")

heat_mat <- heat_df %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.matrix()

n_cols <- 100
heat_pal <- scales::div_gradient_pal(low = "#4575b4", mid = "white", high = "#d73027")(seq(0,1,length.out = n_cols))
heat_breaks <- seq(-1, 1, length.out = n_cols + 1)

# Show in Plots pane
pheatmap(
  heat_mat,
  color = heat_pal,
  breaks = heat_breaks,
  cluster_rows = TRUE, cluster_cols = TRUE,
  border_color = NA, fontsize = 15,
  main = paste0("DL Genes: Gene–Drug Correlations (", corr_method, ")")
)

# Save high-res
pheatmap(
  heat_mat, color = heat_pal, breaks = heat_breaks,
  cluster_rows = TRUE, cluster_cols = TRUE,
  border_color = NA, fontsize = 18,
  main = paste0("DL Genes: Gene–Drug Correlations (", corr_method, ")"),
  filename = "cor_heatmap_TOP_GENES.png", width = png_width/100, height = png_height/100, dpi = png_res
)
pheatmap(
  heat_mat, color = heat_pal, breaks = heat_breaks,
  cluster_rows = TRUE, cluster_cols = TRUE,
  border_color = NA, fontsize = 12,
  main = paste0("DL Genes: Gene–Drug Correlations (", corr_method, ")"),
  filename = "cor_heatmap_TOP_GENES.pdf", width = pdf_width, height = pdf_height
)

dim(heat_mat)

head(heat_mat)

sum(is.na(heat_mat))

sum(is.finite(as.matrix(heat_mat)))

pheatmap(
  heat_mat,
  color = heat_pal,
  breaks = heat_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  
  main = "Top Gene and Drug Correlation",
  
  fontsize = 12,          # row and column labels
  fontsize_row = 12,
  fontsize_col = 12,
  fontsize_number = 12,
  
  fontfamily = "Times New Roman",
  
  angle_col = 45
)
warnings()


# ==========================================================
# Volcano plots only — DL genes validation (Expression ↔ GDSC)
# Single-drug volcano + All-drugs faceted volcano
# ==========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

# ---------------------------
# User parameters (edit here)
# ---------------------------

# Input paths
exp_path_raw    <- "Expression_Transposed_CellNames_RENAMED.csv"
exp_path_clean  <- "Expression_clean_harmonized.csv"  # if present, will be used
gdsc_path       <- "GDSC_AUC_with_cellname_modelID.csv"
model_path      <- "Model.csv"
top_genes_path  <- "TOP GENES.txt"                    # your DL-predicted genes

# Use already harmonized expression file if present?
use_harmonized_expression <- TRUE

# Analyze only your DL genes (TOP GENES.txt)?
use_top_genes_only <- TRUE     # set FALSE to use all genes in expression

# Correlation options
corr_method   <- "spearman"    # "pearson" or "spearman"
min_samples   <- 30            # minimum overlapping cell lines per gene–drug correlation
alpha_fdr     <- 0.05          # FDR threshold for significance highlight
effect_cut    <- 0.25          # |r| threshold for significance highlight

# Volcano (single drug)
volcano_drug  <- NULL          # e.g., "Cisplatin_AUC"; NULL -> auto-select by coverage
label_top_n   <- 15            # label top-N significant genes by |r|

# Volcano (all drugs; stacked/faceted)
facet_ncol    <- 1             # 1 for a tall stack; 2–3 for grid layout
points_alpha  <- 0.75
points_size   <- 1.6

# AUC handling:
# Default: lower AUC = more sensitive (negative r => higher expression -> more sensitivity).
# If you prefer visuals where "higher value = more sensitive", set invert_auc = TRUE.
invert_auc    <- FALSE

# File output sizes
png_width  <- 2200; png_height <- 1800; png_res <- 250  # px
pdf_width  <- 10;   pdf_height <- 8                     # inches

# ---------------------------
# Helpers
# ---------------------------

normalize_name <- function(x) {
  x <- trimws(as.character(x))
  x <- ifelse(is.na(x) | x == "NA", NA_character_, x)
  ifelse(is.na(x), NA_character_, str_replace_all(x, "[^0-9A-Za-z]+", ""))
}

first_non_na <- function(...) {
  vals <- c(...)
  if (all(is.na(vals))) NA_character_ else vals[which(!is.na(vals))[1]]
}

# Robust TOP GENES parser (handles numbering, BOM, blanks)
extract_symbol_safe <- function(line) {
  if (is.na(line)) return(NA_character_)
  line <- trimws(line)
  if (line == "") return(NA_character_)
  line <- sub("^\uFEFF", "", line)                          # strip BOM
  line <- sub("^\\s*\\d+\\s*[:.)-]*\\s*", "", line)         # strip leading numbers/bullets
  toks <- unlist(strsplit(line, "[,\\t\\s]+"))
  toks <- toks[nzchar(toks)]
  if (length(toks) == 0) return(NA_character_)
  sym <- toks[length(toks)]                                 # take last token
  sym <- gsub("[^A-Za-z0-9_\\-]", "", sym)
  if (sym == "") NA_character_ else sym
}

# Efficient correlation of selected genes vs one drug
cor_gene_vs_drug <- function(expr_mat, drug_vec, method = "spearman", min_n = 30) {
  ok_drug <- !is.na(drug_vec)
  if (sum(ok_drug) < min_n) {
    return(tibble(gene = colnames(expr_mat), r = NA_real_, p = NA_real_, n = 0L))
  }
  X <- expr_mat[ok_drug, , drop = FALSE]
  y <- drug_vec[ok_drug]
  
  n_vec <- colSums(!is.na(X) & !is.na(y))
  keep  <- n_vec >= min_n
  res <- tibble(gene = colnames(expr_mat), r = NA_real_, p = NA_real_, n = n_vec)
  if (!any(keep)) return(res)
  
  Xk <- X[, keep, drop = FALSE]
  r_vec <- suppressWarnings(cor(Xk, y, use = "pairwise.complete.obs", method = method))
  nk   <- n_vec[keep]
  r2   <- pmin(pmax(r_vec, -0.999999), 0.999999) # guard ±1
  tval <- r2 * sqrt((nk - 2) / (1 - r2^2))
  pval <- 2 * pt(abs(tval), df = pmax(nk - 2, 1), lower.tail = FALSE)
  
  res$r[keep] <- as.numeric(r_vec)
  res$p[keep] <- as.numeric(pval)
  res
}

# ---------------------------
# Read DL TOP GENES (if used)
# ---------------------------

top_genes <- character(0)
if (use_top_genes_only) {
  if (!file.exists(top_genes_path)) {
    stop("use_top_genes_only=TRUE but 'TOP GENES.txt' not found at: ", top_genes_path)
  }
  top_raw <- readLines(top_genes_path, warn = FALSE, encoding = "UTF-8")
  parsed  <- unlist(lapply(top_raw, extract_symbol_safe), use.names = FALSE)
  top_genes <- unique(na.omit(parsed))
  message("Loaded ", length(top_genes), " DL-predicted genes from: ", top_genes_path)
}

# ---------------------------
# Load & harmonize: Model map
# ---------------------------

model_df <- read_csv(model_path, show_col_types = FALSE)
stopifnot("ModelID" %in% names(model_df))
has_stripped <- "StrippedCellLineName" %in% names(model_df)
has_cellname <- "CellLineName" %in% names(model_df)

model_map <- model_df %>%
  mutate(NormalizedName = dplyr::case_when(
    has_stripped ~ normalize_name(StrippedCellLineName),
    has_cellname ~ normalize_name(CellLineName),
    TRUE ~ NA_character_
  )) %>%
  select(ModelID, NormalizedName) %>%
  distinct()

# ---------------------------
# Load & harmonize: GDSC
# ---------------------------

gdsc_full <- read_csv(gdsc_path, show_col_types = FALSE, guess_max = 1e5) %>%
  mutate(
    cellname_modelID_clean = if ("cellname_modelID" %in% names(.)) normalize_name(cellname_modelID) else NA_character_,
    final_cellname_clean   = if ("final_cellname"   %in% names(.)) normalize_name(final_cellname)   else NA_character_,
    CellLineName_clean     = if ("CellLineName"     %in% names(.)) normalize_name(CellLineName)     else NA_character_
  )

if ("ModelID" %in% names(gdsc_full)) {
  gdsc_full <- gdsc_full %>%
    left_join(model_map, by = "ModelID") %>%
    rename(from_model_map_clean = NormalizedName)
} else {
  gdsc_full <- gdsc_full %>% mutate(from_model_map_clean = NA_character_)
}

priority_cols <- c("cellname_modelID_clean", "final_cellname_clean", "CellLineName_clean", "from_model_map_clean")
present_cols  <- priority_cols[priority_cols %in% names(gdsc_full)]
stopifnot(length(present_cols) > 0)

gdsc_full <- gdsc_full %>%
  mutate(HarmonizedName = pmap_chr(across(all_of(present_cols)), first_non_na))

# Numeric drug columns (exclude name/meta fields)
name_cols <- c("GDSC_CellLineName","ModelID","CellLineName","final_cellname","cellname_modelID",
               "cellname_modelID_clean","final_cellname_clean","CellLineName_clean","from_model_map_clean",
               "HarmonizedName")
drug_cols <- setdiff(names(gdsc_full)[sapply(gdsc_full, is.numeric)], name_cols)
if (length(drug_cols) == 0) stop("No numeric drug columns found in GDSC file.")

# Optional invert AUC
if (invert_auc) {
  gdsc_full[drug_cols] <- lapply(gdsc_full[drug_cols], function(z) ifelse(is.na(z), NA_real_, -z))
}

# ---------------------------
# Load & harmonize: Expression
# ---------------------------

if (use_harmonized_expression && file.exists(exp_path_clean)) {
  exp_df <- read_csv(exp_path_clean, show_col_types = FALSE)
  if (!"CellLineName_clean" %in% names(exp_df) && "CellLineName" %in% names(exp_df)) {
    exp_df <- exp_df %>% mutate(CellLineName_clean = normalize_name(CellLineName))
  }
} else {
  exp_df <- read_csv(exp_path_raw, show_col_types = FALSE)
  stopifnot("CellLineName" %in% names(exp_df))
  exp_df <- exp_df %>% mutate(CellLineName_clean = normalize_name(CellLineName))
}

# De-duplicate cell lines
exp_df <- exp_df %>%
  filter(!is.na(CellLineName_clean)) %>%
  distinct(CellLineName_clean, .keep_all = TRUE)

# ---------------------------
# Join Expression with GDSC by harmonized name
# ---------------------------

joined <- exp_df %>%
  inner_join(
    gdsc_full %>% select(HarmonizedName, all_of(drug_cols)),
    by = c("CellLineName_clean" = "HarmonizedName")
  )

if (nrow(joined) < min_samples) {
  stop("Too few overlapping cell lines after join. Check naming and files.")
}

# ---------------------------
# Build expression matrix (+ optional DL gene filter)
# ---------------------------

expr_name_cols  <- c("CellLineName", "CellLineName_clean")
expr_candidates <- setdiff(names(joined), c(expr_name_cols, drug_cols))
expr_candidates <- expr_candidates[sapply(joined[expr_candidates], is.numeric)]

if (use_top_genes_only) {
  found_genes <- intersect(expr_candidates, top_genes)
  if (length(found_genes) == 0) {
    # case-insensitive fallback
    map_ci <- setNames(expr_candidates, toupper(expr_candidates))
    want   <- unique(toupper(top_genes))
    found_ci <- unique(na.omit(map_ci[intersect(names(map_ci), want)]))
    found_genes <- found_ci
  }
  if (length(found_genes) == 0) {
    stop("None of the TOP GENES were found in expression columns. Check symbols/columns.")
  }
  expr_cols <- found_genes
} else {
  expr_cols <- expr_candidates
}

expr_mat <- as.matrix(joined[expr_cols])
drug_mat <- joined[drug_cols]

# ---------------------------
# Compute correlations for ALL drugs (volcano needs all)
# ---------------------------

all_cor <- vector("list", length(drug_cols))
names(all_cor) <- drug_cols

for (d in drug_cols) {
  dv <- joined[[d]]
  cor_res <- cor_gene_vs_drug(expr_mat, dv, method = corr_method, min_n = min_samples) %>%
    mutate(drug = d)
  all_cor[[d]] <- cor_res
}

cor_long <- bind_rows(all_cor) %>%
  group_by(drug) %>%
  mutate(FDR = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(sig = FDR < alpha_fdr & abs(r) >= effect_cut,
         neglog10FDR = -log10(pmax(FDR, 1e-300)))

# Save table (optional)
write_csv(cor_long, if (use_top_genes_only) "correlation_long_VOLCANO_TOP_GENES.csv" else "correlation_long_VOLCANO_ALL_GENES.csv")

# ---------------------------
# Volcano: single drug
# ---------------------------

if (is.null(volcano_drug)) {
  # auto-pick: drug with most non-NA entries
  non_na_by_drug <- sapply(drug_mat, function(z) sum(!is.na(z)))
  volcano_drug <- names(sort(non_na_by_drug, decreasing = TRUE))[1]
}
stopifnot(volcano_drug %in% cor_long$drug)

vol_df <- cor_long %>%
  filter(drug == volcano_drug, !is.na(r), n >= min_samples)

# Labels: top |r| among significant
lab_genes <- vol_df %>%
  filter(sig) %>%
  slice_max(order_by = abs(r), n = label_top_n, with_ties = FALSE) %>%
  pull(gene)

p_volcano <- ggplot(vol_df, aes(x = r, y = neglog10FDR)) +
  geom_point(aes(color = sig), alpha = 0.85, size = 1.9) +
  scale_color_manual(values = c(`TRUE` = "#d73027", `FALSE` = "grey60")) +
  geom_vline(xintercept = c(-effect_cut, effect_cut), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(alpha_fdr), linetype = 2, color = "grey50") +
  ggrepel::geom_text_repel(
    data = subset(vol_df, gene %in% lab_genes),
    aes(label = gene), size = 3, max.overlaps = Inf, min.segment.length = 0
  ) +
  labs(
    title = paste0("Volcano: ", volcano_drug, " (", corr_method, ")"),
    subtitle = if (invert_auc) "AUC inverted (higher value = more sensitive)"
    else "Lower AUC = more sensitive",
    x = "Correlation (r)",
    y = expression(-log10),
    color = "Significant"
  ) +
  theme_minimal(base_size = 12)

# Show in Plots pane
print(p_volcano)

# Save
ggsave(paste0("VOLCANO_", make.names(volcano_drug),
              if (use_top_genes_only) "_TOP_GENES.png" else "_ALL_GENES.png"),
       p_volcano, width = png_width / png_res, height = png_height / png_res, dpi = png_res, units = "in")
ggsave(paste0("VOLCANO_", make.names(volcano_drug),
              if (use_top_genes_only) "_TOP_GENES.pdf" else "_ALL_GENES.pdf"),
       p_volcano, width = pdf_width, height = pdf_height)

# ---------------------------
# Volcano: ALL drugs (stacked/faceted)
# ---------------------------

# Order panels by coverage
drug_order <- cor_long %>%
  group_by(drug) %>% summarise(m = sum(!is.na(r)), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(drug)

volc_df_all <- cor_long %>%
  mutate(drug = factor(drug, levels = drug_order))

p_volc_stack <- ggplot(volc_df_all, aes(x = r, y = neglog10FDR)) +
  geom_point(aes(color = sig), size = points_size, alpha = points_alpha) +
  scale_color_manual(values = c(`TRUE` = "#d73027", `FALSE` = "grey60")) +
  geom_vline(xintercept = c(-effect_cut, effect_cut), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(alpha_fdr), linetype = 2, color = "grey50") +
  labs(
    title = "Volcano plots for all drugs",
    subtitle = paste0("Spearman r vs –log10(FDR); thresholds: |r|≥", effect_cut, ", FDR<", alpha_fdr,
                      ifelse(invert_auc, " | AUC inverted (higher = more sensitive)", " | lower AUC = more sensitive")),
    x = "Correlation (r)",
    y = expression(-log10),
    color = "Significant"
  ) +
  theme_minimal(base_size = 11) +
  facet_wrap(~ drug, ncol = facet_ncol, scales = "free_y") +
  theme(legend.position = "top")

# Ensure finite r and neglog10FDR only (avoid warnings about dropped rows)
volc_df_all <- cor_long %>%
  filter(!is.na(r), is.finite(r),
         !is.na(neglog10FDR), is.finite(neglog10FDR),
         n >= min_samples) %>%
  mutate(sig = FDR < alpha_fdr & abs(r) >= effect_cut)

# (re-)order panels by coverage (optional)
drug_order <- volc_df_all %>%
  group_by(drug) %>% summarise(m = sum(!is.na(r)), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(drug)
volc_df_all <- volc_df_all %>%
  mutate(drug = factor(drug, levels = drug_order))

print(p_volc_stack)

# ==========================================
# SAFE SAVE for stacked/faceted volcano plots
# - Finite filtering (avoid "Removed rows")
# - Auto pagination (pages of r x c facets)
# - Auto column choice for single-page export
# - Vector PDF (no pixel limit) + PNG pages
# ==========================================

# ==========================================
# SAFE SAVE for stacked/faceted volcano plots
# - Finite filtering (avoid "Removed rows")
# - Auto pagination (pages of r x c facets)
# - Auto column choice for single-page export
# - Vector PDF (no pixel limit) + PNG pages
# ==========================================


unlink("C:/Users/22994327/AppData/Local/Programs/R/R-4.4.2/library/00LOCK", recursive = TRUE, force = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  install.packages("ggforce")
  library(ggforce)   # for facet_wrap_paginate (install.packages("ggforce") if needed)
})

# --- 1) Ensure finite values only (avoid warnings) ---
volc_df_all <- cor_long %>%
  filter(!is.na(r), is.finite(r),
         !is.na(neglog10FDR), is.finite(neglog10FDR),
         n >= min_samples) %>%
  mutate(sig = FDR < alpha_fdr & abs(r) >= effect_cut)

# --- 2) Order panels by coverage (optional) ---
drug_order <- volc_df_all %>%
  group_by(drug) %>% summarise(m = sum(!is.na(r)), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(drug)

volc_df_all <- volc_df_all %>%
  mutate(drug = factor(drug, levels = drug_order))

# --- 3) Core plot (no facets yet) ---
p_base <- ggplot(volc_df_all, aes(x = r, y = neglog10FDR)) +
  geom_point(aes(color = sig), size = points_size, alpha = points_alpha) +
  scale_color_manual(values = c(`TRUE` = "#d73027", `FALSE` = "grey60")) +
  geom_vline(xintercept = c(-effect_cut, effect_cut), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(alpha_fdr), linetype = 2, color = "grey50") +
  labs(
    title = "Volcano plots for all drugs",
    subtitle = paste0("Spearman r vs –log10(FDR); thresholds: |r|≥", effect_cut,
                      ", FDR<", alpha_fdr,
                      ifelse(invert_auc, " | AUC inverted (higher = more sensitive)",
                             " | lower AUC = more sensitive)")),
    x = "Correlation (r)",
    y = expression(-log10),
    color = "Significant"
  ) +
  theme_minimal(base_size = 11)
# If you want comparable x-axis across drugs, uncomment:
# + coord_cartesian(xlim = c(-1, 1))

# --- 4) Multi-page export helper (pagination) ---
save_paginated_volcanos <- function(df, plot_base,
                                    page_nrow = 4, page_ncol = 3,
                                    file_stub = "VOLCANO_ALL_DRUGS",
                                    width_in = pdf_width, height_in = pdf_height,
                                    dpi = png_res) {
  n_panels <- length(levels(df$drug))
  n_pages  <- ceiling(n_panels / (page_nrow * page_ncol))
  message(sprintf("Paginating %d panels into %d page(s) of %dx%d facets.",
                  n_panels, n_pages, page_nrow, page_ncol))
  
  for (i in seq_len(n_pages)) {
    p_page <- plot_base +
      ggforce::facet_wrap_paginate(~ drug, nrow = page_nrow, ncol = page_ncol,
                                   scales = "free_y", page = i) +
      theme(legend.position = "top")
    
    print(p_page)  # Plots pane preview
    
    # Vector PDF page (no pixel limit)
    ggsave(sprintf("%s_p%02d.pdf", file_stub, i), p_page,
           width = width_in, height = height_in)
    
    # PNG page (bounded size)
    ggsave(sprintf("%s_p%02d.png", file_stub, i), p_page,
           width = width_in, height = height_in, units = "in", dpi = dpi)
  }
}

# --- 5) One-page export helper (auto ncol to keep height & pixels under limits) ---
save_singlepage_volcanos <- function(df, plot_base,
                                     min_row_height_in = 0.55,  # vertical density per row
                                     base_margins_in   = 1.0,   # extra top/bottom space
                                     max_in            = 49.5,  # < 50 in guard
                                     max_px            = 49990, # < 50k px device guard
                                     width_in          = pdf_width,
                                     dpi               = png_res,
                                     file_stub         = "VOLCANO_ALL_DRUGS_stack") {
  n_panels <- length(levels(df$drug))
  
  # Try columns from 1 up to sqrt(n_panels) to find a layout that fits both inch & pixel caps
  choose_layout <- function() {
    best <- NULL
    for (nc in 1:max(1, ceiling(sqrt(n_panels)))) {
      n_rows <- ceiling(n_panels / nc)
      h_in   <- max(pdf_height, n_rows * min_row_height_in + base_margins_in)
      h_in   <- min(h_in, max_in)   # cap to <50 inches
      # Pixel check for PNG
      h_px   <- h_in * dpi
      w_px   <- width_in * dpi
      fits_px <- (h_px <= max_px && w_px <= max_px)
      # Accept as soon as we find one that fits
      if (h_in <= max_in && fits_px) return(list(ncol = nc, height_in = h_in))
      # Keep the best (smallest height) even if PNG doesn't fit, for PDF
      if (is.null(best) || h_in < best$height_in) best <- list(ncol = nc, height_in = h_in)
    }
    best
  }
  
  choice <- choose_layout()
  facet_ncol_auto <- choice$ncol
  height_in <- choice$height_in
  
  message(sprintf("Single-page layout → ncol=%d, rows≈%d, height=%.2f in",
                  facet_ncol_auto, ceiling(n_panels / facet_ncol_auto), height_in))
  
  p_faceted <- plot_base +
    facet_wrap(~ drug, ncol = facet_ncol_auto, scales = "free_y") +
    theme(legend.position = "top")
  
  # Save vector PDF (always safe)
  ggsave(paste0(file_stub, ".pdf"), p_faceted,
         width = width_in, height = height_in)
  
  # Save PNG if within pixel bounds; otherwise fall back to pagination
  if ((height_in * dpi) <= max_px && (width_in * dpi) <= max_px) {
    ggsave(paste0(file_stub, ".png"), p_faceted,
           width = width_in, height = height_in, units = "in", dpi = dpi)
  } else {
    message("PNG would exceed device pixel limits → switching to paginated PNGs.")
    save_paginated_volcanos(df, plot_base,
                            page_nrow = 4, page_ncol = 3,
                            file_stub = "VOLCANO_ALL_DRUGS",
                            width_in = width_in, height_in = pdf_height, dpi = dpi)
  }
  
  invisible(p_faceted)
}

# --- 6) Choose which strategy you prefer ---

# (A) Multi-page output (recommended for many drugs)
save_paginated_volcanos(
  df        = volc_df_all,
  plot_base = p_base,
  page_nrow = 4, page_ncol = 3,               # 12 facets per page (tweak as you like)
  file_stub = "VOLCANO_ALL_DRUGS",
  width_in  = pdf_width, height_in = pdf_height, dpi = png_res
)

# (B) Single-page attempt (auto-adjust columns & height; falls back to pagination for PNG)
save_singlepage_volcanos(
  df         = volc_df_all,
  plot_base  = p_base,
  min_row_height_in = 0.55,
  base_margins_in   = 1.0,
  max_in     = 49.5,                            # keep < 50 inches
  max_px     = 49990,                           # keep < 50,000 px
  width_in   = pdf_width,
  dpi        = png_res,
  file_stub  = "VOLCANO_ALL_DRUGS_stack"
)

# Tip: If you only want one of (A) or (B), comment the other out.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
})

# --- 1) Identify the drugs shown on your heatmap ---
# If you have heat_df (wide: gene + drug columns):
if (exists("heat_df")) {
  heatmap_drugs <- setdiff(colnames(heat_df), "gene")
} else if (exists("heat_mat")) {
  heatmap_drugs <- colnames(heat_mat)
} else {
  stop("Neither heat_df nor heat_mat found. Run the heatmap selection first.")
}

if (length(heatmap_drugs) == 0) stop("No drugs detected in the heatmap object.")

message("Volcano will be plotted for heatmap drugs: ",
        paste(heatmap_drugs, collapse = ", "))

# --- 2) Filter cor_long to heatmap drugs & compute plotting fields ---
volc_heat <- cor_long %>%
  filter(drug %in% heatmap_drugs, !is.na(r), !is.na(FDR), n >= min_samples) %>%
  mutate(
    sig         = (FDR < alpha_fdr & abs(r) >= effect_cut),
    neglog10FDR = -log10(pmax(FDR, 1e-300))
  )

# Optional: order facets by number of significant hits (or by coverage)
drug_order <- volc_heat %>%
  group_by(drug) %>%
  summarise(sig_n = sum(sig, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(sig_n)) %>%
  pull(drug)

volc_heat <- volc_heat %>%
  mutate(drug = factor(drug, levels = drug_order))

# --- 3) Faceted volcano (heatmap drugs only) ---
label_top_n <- 10  # change to 0 if you don't want labels in the faceted figure

lab_df <- volc_heat %>%
  filter(sig) %>%
  group_by(drug) %>%
  slice_max(order_by = abs(r), n = label_top_n, with_ties = FALSE) %>%
  ungroup()

p_volc_heat <- ggplot(volc_heat, aes(x = r, y = neglog10FDR)) +
  geom_point(aes(color = sig), alpha = 0.8, size = 1.7) +
  scale_color_manual(values = c(`TRUE` = "#d73027", `FALSE` = "grey60")) +
  geom_vline(xintercept = c(-effect_cut, effect_cut), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(alpha_fdr), linetype = 2, color = "grey50") +
  ggrepel::geom_text_repel(
    data = lab_df, aes(label = gene),
    size = 3, max.overlaps = Inf, min.segment.length = 0, box.padding = 0.3
  ) +
  labs(
    title = "Volcano plots for heatmap drugs only",
    subtitle = paste0("Spearman r vs –log10(FDR); thresholds: |r|≥", effect_cut,
                      ", FDR<", alpha_fdr,
                      ifelse(exists("invert_auc") && invert_auc,
                             " (AUC inverted: higher = more sensitive)",
                             " (lower AUC = more sensitive)")),
    x = "Correlation (r)",
    y = expression(-log10),
    color = "Significant"
  ) +
  theme_minimal(base_size = 11) +
  facet_wrap(~ drug, ncol = 3, scales = "free_y")  # adjust columns to taste

# Show in Plots pane
print(p_volc_heat)

# --- 4) Save a standard-sized figure (journal-friendly) ---
# Adjust width/height depending on number of panels
n_panels <- length(levels(volc_heat$drug))
ncol_fac <- 3
nrow_fac <- ceiling(n_panels / ncol_fac)
fig_w <- 10
fig_h <- max(6, 2.6 * nrow_fac)   # ~2.6 inches per row

ggsave("VOLCANO_heatmap_drugs_facet.pdf", p_volc_heat, width = fig_w, height = fig_h)
ggsave("VOLCANO_heatmap_drugs_facet.png", p_volc_heat, width = fig_w, height = fig_h, dpi = 300)