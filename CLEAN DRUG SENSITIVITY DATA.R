setwd(setwd("C:/DESKTOP/BIOINFORMATICS VALIDATION"))

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# ---- 1) Read both files (semicolon-delimited) ----
f1 <- "PANCANCER_IC_Wed Mar  4 10_31_25 2026_DRUG SENSITIVITY 1.csv" # (GDSC2)
f2 <- "PANCANCER_IC_Wed Mar  4 10_34_19 2026_DRUG SENSITIVITY 2.csv" # (GDSC1)

g1 <- read_delim(f1, delim = ";", col_types = cols(.default = "c"))
g2 <- read_delim(f2, delim = ";", col_types = cols(.default = "c"))

# ---- 2) Standardize column names (trim spaces; collapse multiple spaces) ----
clean_names <- function(x) {
  x %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}
names(g1) <- clean_names(names(g1))
names(g2) <- clean_names(names(g2))

# ---- 3) Verify required columns exist (Drug Name, Cell Line Name, AUC) ----
needed <- c("Drug Name","Cell Line Name","AUC")
stopifnot(all(needed %in% names(g1)) | all(needed %in% names(g2)))

# ---- 4) Select only what we need and coerce AUC to numeric ----
pick_cols <- function(df) {
  keep <- intersect(c("Drug Name","Cell Line Name","AUC","Dataset Version"), names(df))
  df %>%
    select(all_of(keep)) %>%
    mutate(
      `Drug Name`      = as.character(`Drug Name`),
      `Cell Line Name` = as.character(`Cell Line Name`),
      AUC              = suppressWarnings(as.numeric(AUC))
    )
}

g1_s <- pick_cols(g1)
g2_s <- pick_cols(g2)

# ---- 5) Combine both files ----
gdsc_long <- bind_rows(g1_s, g2_s) %>%
  distinct(`Drug Name`, `Cell Line Name`, AUC, .keep_all = TRUE)

# Optional: If you only want GDSC1 OR GDSC2, uncomment one of these:
# gdsc_long <- gdsc_long %>% filter(`Dataset Version` == "GDSC1")
# gdsc_long <- gdsc_long %>% filter(`Dataset Version` == "GDSC2")

# ---- 6) Resolve duplicate (Cell line, Drug) by averaging AUC ----
gdsc_long_agg <- gdsc_long %>%
  group_by(`Cell Line Name`, `Drug Name`) %>%
  summarise(AUC = mean(AUC, na.rm = TRUE), .groups = "drop")

write_csv(gdsc_long_agg, "GDSC_all_AUC_matrix.csv")

#Transpose cell lines to rows and drugs to columns

library(tidyverse)

# Read without altering column names
df <- read.csv("GDSC_all_AUC_matrix.csv", check.names = FALSE)

# Sanity check: see the exact names as in the file
print(names(df))

# Pivot: rows = Cell Line Name, columns = Drug Name, values = AUC
wide_df <- df %>%
  pivot_wider(
    id_cols    = `Cell Line Name`,
    names_from = `Drug Name`,
    values_from = AUC,
    values_fn  = mean
  )

write.csv(wide_df, "GDSC_AUC_pivot.csv", row.names = FALSE)
cat("Done. Dimensions:", nrow(wide_df), "rows ×", ncol(wide_df), "columns\n")

# Optional: Print size of resulting matrix
cat("Pivoted matrix dimensions:", nrow(wide_df), "rows ×", ncol(wide_df), "columns\n")
