setwd(setwd("C:/DESKTOP/BIOINFORMATICS VALIDATION"))
###############################################
## CLEAN CCLE TPM LOG(TPM+1)
## KEEP ACH MODEL IDs AS TEXT (never numeric)
## Rows = Cell lines; Columns = Genes
###############################################

library(tidyverse)
library(data.table)
library(stringr)

#----------------------------------------------------------
# 1. LOAD DATASET SAFELY (ACH IDs remain character)
#----------------------------------------------------------



library(data.table)
library(tidyverse)

# 1) Load robustly (handle broken quotes/rows). Do NOT set colClasses yet.
expr <- fread(
  "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  header = TRUE,
  fill   = TRUE,   # tolerate uneven rows
  quote  = "",     # ignore broken quotes
  sep    = ",",
  check.names = FALSE
)

# 2) Ensure the first column (Model ID) stays as character.
#    We don't assume the name; we detect it dynamically.
first_col <- names(expr)[1]
expr[[first_col]] <- as.character(expr[[first_col]])

# 3) Make duplicate headers unique to prevent tibble errors later.
names(expr) <- make.unique(names(expr))



# Repair duplicated column names from broken header
names(expr) <- make.unique(names(expr))

# Check first 20 column names
head(names(expr), 20)


#----------------------------------------------------------
# 2. Detect gene columns (they look like "TSPAN6 (7105)")
#----------------------------------------------------------

is_gene <- grepl("\\([0-9]+\\)$", names(expr))
sum(is_gene)  # how many genes detected


#----------------------------------------------------------
# 3. Identify Model ID column
#----------------------------------------------------------

model_id_col <- names(expr)[1]
cat("Model ID column detected:", model_id_col, "\n")


#----------------------------------------------------------
# 4. Clean gene names: "SYMBOL (1234)" → "SYMBOL"
#----------------------------------------------------------

gene_cols <- names(expr)[is_gene]

gene_map <- tibble(
  raw = gene_cols,
  SYMBOL = str_remove(gene_cols, "\\s*\\(.*\\)$"),
  ENTREZ = str_extract(gene_cols, "(?<=\\()[0-9]+(?=\\))")
)

# Replace raw names with cleaned SYMBOL names
names(expr)[match(gene_cols, names(expr))] <- gene_map$SYMBOL


#----------------------------------------------------------
# 5. Convert gene expression columns to numeric
#----------------------------------------------------------

expr <- expr %>%
  mutate(across(all_of(gene_map$SYMBOL),
                ~ suppressWarnings(as.numeric(.x))))


#----------------------------------------------------------
# 6. Reduce to MODEL ID + GENE columns only
#----------------------------------------------------------

expr_small <- expr %>%
  select(all_of(model_id_col), all_of(gene_map$SYMBOL))


#----------------------------------------------------------
# 7. TRANSPOSE → Rows = GENES, Columns = ACH Model IDs
#----------------------------------------------------------

expr_matrix <- expr_small %>%
  column_to_rownames(var = model_id_col) %>%
  as.matrix()

expr_genes_as_rows <- t(expr_matrix)

expr_final <- expr_genes_as_rows %>%
  as_tibble(rownames = "Gene")


#----------------------------------------------------------
# 8. Inspect and Save final output
#----------------------------------------------------------

dim(expr_final)
head(expr_final[, 1:10])

write_csv(expr_final,
          "OmicsExpressionTPMLogp1HumanProteinCodingGenes_TRANSPOSED.csv")

message("DONE: Transposed matrix saved. ACH model IDs preserved as text.")



cols_to_remove <- c(
  "SequencingID",
  "IsDefaultEntryForModel",
  "ModelConditionID",
  "IsDefaultEntryForMC"
)

expr_clean <- expr_clean %>%
  dplyr::select(-any_of(cols_to_remove))

# Check remaining columns
names(expr_clean)

write_csv(expr_clean, "Expression_Clean.csv")



library(dplyr)

model_id_col <- names(expr_clean)[1]
expr_clean[[model_id_col]] <- as.character(expr_clean[[model_id_col]])

# 1) Identify gene columns (everything except ModelID)
gene_cols <- setdiff(names(expr_clean), model_id_col)

# 2) Ensure gene columns are numeric
expr_clean <- expr_clean %>% mutate(across(all_of(gene_cols), ~ suppressWarnings(as.numeric(.x))))

# 3) Aggregate by ACH (mean across duplicates, na.rm = TRUE)
expr_clean <- expr_clean %>%
  group_by(.data[[model_id_col]]) %>%
  summarise(across(all_of(gene_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 4) Now transpose (genes = rows, ACH = columns)
expr_t <- expr_clean %>%
  tibble::column_to_rownames(var = model_id_col) %>%
  as.matrix() %>%
  t() %>%
  tibble::as_tibble(rownames = "Gene")


write_csv(expr_clean, "expr_t.csv")

# 5) Concert ModelID to CellName
expr_t <- read.csv("Expr_t.csv", check.names = FALSE)

model <- read.csv("Model.csv", check.names = FALSE)
head(model[, c("ModelID", "CellLineName")])
ach_ids <- colnames(expr_t)[-1]   # skip the Gene column
id_map <- model %>%
  dplyr::filter(ModelID %in% ach_ids) %>%
  dplyr::select(ModelID, CellLineName)

new_names <- colnames(expr_t)

for (i in 2:length(new_names)) {
  ach <- new_names[i]
  cell <- id_map$CellLineName[id_map$ModelID == ach]
  if (length(cell) == 1 && !is.na(cell)) {
    new_names[i] <- paste0(cell, "_", ach)
  }
}

colnames(expr_t) <- new_names

write.csv(expr_t, "Expression_Transposed_CellNames.csv", row.names = FALSE)



