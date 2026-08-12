###############################################################
## FINAL CLEAN SCRIPT
## Drug Sensitivity → Differential Expression → GSEA → Panels
###############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
  library(fgsea)
  library(patchwork)
  library(scales)
  library(stringr)
})

###############################################################
# 1. Helper functions
###############################################################

# Canonical ID
canon <- function(x){
  x %>% toupper() %>% str_replace_all("[^A-Z0-9]", "") %>% str_trim()
}

# Welch t-test per gene
row_t_welch <- function(mat, grp1, grp2){
  if(length(grp1)<2 || length(grp2)<2) return(NULL)
  x <- mat[grp1,,drop=FALSE]
  y <- mat[grp2,,drop=FALSE]
  
  m1 <- colMeans(x, na.rm=TRUE)
  m2 <- colMeans(y, na.rm=TRUE)
  v1 <- apply(x,2,var,na.rm=TRUE)
  v2 <- apply(y,2,var,na.rm=TRUE)
  n1 <- colSums(!is.na(x))
  n2 <- colSums(!is.na(y))
  se <- sqrt(v1/n1 + v2/n2)
  t  <- (m1 - m2) / se
  df <- (v1/n1 + v2/n2)^2 / ((v1^2/((n1^2)*(n1-1))) + (v2^2/((n2^2)*(n2-1))))
  p  <- 2 * pt(-abs(t), df)
  
  tibble(gene = names(m1), t=t, p=p, logFC=m1-m2)
}

# Ranking metric for GSEA
rank_genes <- function(t_stat, pval){
  pval[pval==0] <- 1e-300
  score <- t_stat * -log10(pval)
  score[is.na(score)] <- 0
  sort(score, decreasing = TRUE)
}

# Panel plotting (NES strip + bubble)
make_panel <- function(df, show_top=12){
  
  df_sel <- df %>% arrange(padj, desc(abs(NES))) %>% head(show_top)
  if(nrow(df_sel)==0) return(ggplot()+theme_void())
  
  df_sel <- df_sel %>%
    mutate(label = factor(label, levels = rev(label)),
           sig_cat = if_else(padj<=0.05,"padj ≤ 0.05","padj > 0.05"))
  
  col_strip <- c("padj > 0.05"="#ff6b81","padj ≤ 0.05"="#12d3a6")
  
  p_left <- ggplot(df_sel, aes(NES, label, fill=sig_cat)) +
    geom_col() +
    geom_point(
      data = df_sel %>% filter(padj<=0.05),
      aes(x = NES + 0.05*sign(NES), y = label),
      shape=8, size=2.5, color="black"
    ) +
    scale_fill_manual(values=col_strip) +
    ggtitle(unique(df_sel$drug)) +
    theme_minimal(base_size=10) +
    theme(legend.position="top")
  
  p_right <- ggplot(df_sel, aes(ratio, label)) +
    geom_point(aes(size=k, color=padj), alpha=0.9) +
    scale_color_gradientn(colours=c("blue","red","orange"),
                          trans="log10", name="padj") +
    scale_size(range=c(3,8), name="k") +
    theme_minimal(base_size=10) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  p_left + p_right + plot_layout(widths=c(1,2))
}

###############################################################
# 2. Load & Clean Expression (canonicalize + collapse duplicates)
###############################################################

expr_raw <- read_csv("Expression_clean_harmonized.csv", show_col_types = FALSE)
id_cols <- c("CellLineName_clean","final_cellname","CellLineName")
id_col  <- intersect(id_cols, names(expr_raw))[1]

expr_tbl <- expr_raw %>%
  rename(Cell = !!id_col) %>%
  mutate(CellCanon = canon(Cell))

expr_collapsed <- expr_tbl %>%
  select(CellCanon, where(is.numeric)) %>%
  group_by(CellCanon) %>%
  summarise(across(everything(), ~mean(.x,na.rm=TRUE)), .groups="drop")

expr_df <- as.data.frame(expr_collapsed)
rownames(expr_df) <- expr_df$CellCanon
expr_df$CellCanon <- NULL

expr_df[!is.finite(as.matrix(expr_df))] <- NA_real_
expr_df <- expr_df[, colSums(is.na(expr_df)) < nrow(expr_df), drop=FALSE]

expr_log <- log1p(expr_df)

cat("Expression samples after cleaning:", nrow(expr_log), "\n")

###############################################################
# 3. Load & Clean AUC (canonicalize + collapse duplicates)
###############################################################

auc_raw <- read_csv("GDSC_AUC_with_cellname_modelID.csv", show_col_types = FALSE)
auc_id_cols <- c("final_cellname","CellLineName","CellLineName_clean")
auc_id <- intersect(auc_id_cols, names(auc_raw))[1]

auc_tbl <- auc_raw %>%
  rename(Cell = !!auc_id) %>%
  mutate(CellCanon = canon(Cell))

drugs_12 <- c("AZD8186","Avagacestat","BMS-536924","Afatinib","Foretinib",
              "AZD4547","AZD7762","AZD6738","Gemcitabine","Docetaxel",
              "Camptothecin","5-Fluorouracil")

available_drugs <- intersect(drugs_12, names(auc_tbl))

auc_c <- auc_tbl %>%
  select(CellCanon, all_of(available_drugs)) %>%
  group_by(CellCanon) %>%
  summarise(across(everything(), ~mean(.x,na.rm=TRUE)), .groups="drop")

cat("AUC samples after cleaning:", nrow(auc_c), "\n")

###############################################################
# 4. Sensitive / Resistant grouping (AUC tertiles)
###############################################################

get_groups <- function(drug){
  v <- auc_c[[drug]]
  cells <- auc_c$CellCanon
  ok <- !is.na(v)
  v <- v[ok]; cells <- cells[ok]
  qs <- quantile(v, c(0.33,0.67), na.rm=TRUE)
  sens <- cells[v <= qs[1]]
  resi <- cells[v >= qs[2]]
  
  sens <- intersect(sens, rownames(expr_log))
  resi <- intersect(resi, rownames(expr_log))
  
  list(sens=sens, resi=resi)
}

###############################################################
# 5. Compute genome-wide t-statistics per drug
###############################################################

ranks_list <- list()
for(drug in available_drugs){
  
  g <- get_groups(drug)
  if(length(g$sens)<2 || length(g$resi)<2){
    message("Skipping ", drug, " (insufficient samples).")
    next
  }
  
  res <- row_t_welch(expr_log, g$sens, g$resi)
  res <- res %>% filter(is.finite(t), is.finite(p))
  
  ranks <- rank_genes(res$t, res$p)
  names(ranks) <- res$gene
  
  ranks_list[[drug]] <- ranks
}

cat("Drugs with valid rank lists:", names(ranks_list), "\n")

stopifnot(length(ranks_list) > 0)

###############################################################
# 6. GSEA (Hallmark)
###############################################################

msig_h <- msigdbr(species="Homo sapiens", category="H")
hallmark_list <- split(msig_h$gene_symbol, msig_h$gs_name)

run_gsea <- function(drug, ranks){
  fg <- fgsea(
    pathways = hallmark_list,
    stats = ranks,
    nperm = 10000,
    minSize = 5,
    maxSize = 5000
  )
  
  k <- vapply(fg$leadingEdge, length, numeric(1))
  label <- gsub("^HALLMARK_","",fg$pathway)
  label <- gsub("_"," ",label) %>% str_to_title()
  
  tibble(
    drug = drug,
    pathway = fg$pathway,
    label = label,
    NES = fg$NES,
    padj = fg$padj,
    size = fg$size,
    k = k,
    ratio = k/fg$size,
    direction = if_else(fg$NES>=0,"Up (↑)","Down (↓)")
  )
}

gsea_tbl <- map2_df(names(ranks_list), ranks_list, run_gsea)

###############################################################
# 7. Make multi-panel figure (3 × 4 grid)
###############################################################

# Keep your preferred order; restrict to drugs present in gsea_tbl
grid_order <- intersect(
  c("AZD8186","Avagacestat","BMS-536924","Afatinib",
    "Foretinib","AZD4547","AZD7762","AZD6738",
    "Gemcitabine","Docetaxel","Camptothecin","5-Fluorouracil"),
  unique(gsea_tbl$drug)
)

# Build one panel per drug in that order
panels <- lapply(grid_order, function(d) {
  make_panel(dplyr::filter(gsea_tbl, drug == d), show_top = 12)
})

# If fewer than 12 panels produced (e.g., a drug lacked enough S/R samples),
# pad with empty placeholders to keep a strict 3 × 4 layout.
if (length(panels) < 12) {
  panels <- c(panels, rep(list(ggplot() + theme_void()), 12 - length(panels)))
}

# Assemble as 3 columns × 4 rows
full_plot_3x4 <-
  patchwork::wrap_plots(panels, ncol = 3, nrow = 4, byrow = TRUE) +
  patchwork::plot_annotation(title = "GSEA — Hallmark Pathways (3 × 4)") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save (larger page to fit 12 panels clearly)
dir.create("GSEA_panels", showWarnings = FALSE)

ggsave("GSEA_panels/Hallmark_panels_3x4.pdf",
       full_plot_3x4, width = 16, height = 12, units = "in")

ggsave("GSEA_panels/Hallmark_panels_3x4.png",
       full_plot_3x4, width = 16, height = 12, units = "in",
       dpi = 300, bg = "white")