#!/usr/bin/env Rscript

# Merge Fisher antigen-level results (lev0/1/2) and keep only significant rows
# Comments in English only

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# --- INPUT BASE DIRECTORY (clean version) ---
base_dir <- "/doctorai/chiarba/AbAg_database/clean"

# --- INPUT FILES (from step 3 / Python Fisher) ---
files <- c(
  lev0 = file.path(base_dir, "fisher_results_annot1_antigen_lev0.csv"),
  lev1 = file.path(base_dir, "fisher_results_annot1_antigen_lev1.csv"),
  lev2 = file.path(base_dir, "fisher_results_annot1_antigen_lev2.csv")
)

# --- READ EACH LEVEL AND ADD 'lev' COLUMN ---
dfs <- lapply(names(files), function(nm) {
  df <- read_csv(files[[nm]], show_col_types = FALSE, guess_max = 100000)
  df$lev <- nm
  df
})

# Row-bind all levels
all_levels <- bind_rows(dfs)

# Sanity check for p_adj
if (!"p_adj" %in% names(all_levels)) {
  stop("Column 'p_adj' not found. Please check your CSV headers.")
}

# Ensure p_adj is numeric
all_levels <- all_levels %>%
  mutate(p_adj = as.numeric(p_adj))

# --- FILTER SIGNIFICANT RESULTS AND SORT ---
sig <- all_levels %>%
  filter(!is.na(p_adj) & p_adj < 0.05) %>%
  arrange(lev, p_adj)

# --- OUTPUT (clean, non-overwriting name) ---
out_path <- file.path(base_dir, "fisher_results_annot1_antigen_all_levels_sig.csv")
write_csv(sig, out_path)

# --- SUMMARY ON STDOUT ---
cat(
  "Merged rows:", nrow(all_levels), "\n",
  "Significant rows (p_adj < 0.05):", nrow(sig), "\n",
  "Saved to:", out_path, "\n"
)
