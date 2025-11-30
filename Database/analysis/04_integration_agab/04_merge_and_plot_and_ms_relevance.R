#!/usr/bin/env Rscript

# Merge Fisher antigen-level results (lev0/1/2) and keep only significant rows
# Comments in English only

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# --- INPUT BASE DIRECTORY (clean version) ---
base_dir <- "/doctorai/chiarba/AbAg_database/aggregating_test"

# --- INPUT FILES (from step 3 / Python Fisher) ---
files <- c(
  lev0 = file.path(base_dir, "fisher_results_annot2_antigen_lev0.csv"),
  lev1 = file.path(base_dir, "fisher_results_annot2_antigen_lev1.csv"),
  lev2 = file.path(base_dir, "fisher_results_annot2_antigen_lev2.csv")
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
out_path <- file.path(base_dir, "fisher_results_annot2_antigen_all_levels_sig.csv")
write_csv(sig, out_path)

# --- SUMMARY ON STDOUT ---
cat(
  "Merged rows:", nrow(all_levels), "\n",
  "Significant rows (p_adj < 0.05):", nrow(sig), "\n",
  "Saved to:", out_path, "\n"
)

### add column sig for MS and plotting

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

sig_ms <- fread("/doctorai/chiarba/AbAg_database/aggregating_test/fisher_results_annot2_antigen_all_levels_sig.csv")

# Ensure we have a lowercased version of group
sig_ms <- sig_ms %>%
  mutate(group_chr = as.character(group),
         group_lc  = tolower(group_chr))

# -----------------------------
# MS-relevant antigens tag
# -----------------------------

# Name-based lists for EBV / CMV / Akkermansia / human brain

ebv_ms_names <- c(
  "Epstein-Barr nuclear antigen 1",
  "Epstein-Barr nuclear antigen 2",
  "Epstein-Barr nuclear antigen 6",
  "Epstein-Barr nuclear antigen leader protein",
  "Envelope glycoprotein GP350",
  "G-protein coupled receptor BULF1"
)

cmv_ms_names <- c(
  "Immediate early protein IE1",
  "Serine/threonine protein kinase UL97"
)

akk_ms_names <- c(
  "Amuc_1100"
)


human_brain_names <- c(
  "Neural cell adhesion molecule L1",
  "Neural cell adhesion molecule 1",
  "Microtubule-associated protein tau",
  "Myelin-oligodendrocyte glycoprotein",
  "Sphingomyelinase",
  "Acid sphingomyelinase-like phosphodiesterase 3b"
)

# EBV / CMV / Akkermansia regex (pattern-based)
rx_ebv_ms <- regex(
  "epstein[- ]?barr|\\bebv\\b|gp350|balf1|bulf1",
  ignore_case = TRUE
)

rx_cmv_ms <- regex(
  "cytomegalovirus|\\bhcmv\\b|\\bcmv\\b|immediate early protein ie1|ul97",
  ignore_case = TRUE
)

rx_akk_ms <- regex(
  "akkermansia|amuc_?1100|\\bamuc\\b",
  ignore_case = TRUE
)

# Human brain-related antigens (heuristic list + explicit names above)
brain_terms <- c(
  "myelin basic protein",
  "myelin oligodendrocyte glycoprotein",
  "\\bmog\\b",
  "proteolipid protein 1",
  "\\bplp1\\b",
  "neurofascin",
  "contactin",
  "neurofilament",
  "glial fibrillary acidic protein",
  "\\bgfap\\b",
  "synuclein",
  "huntingtin",
  "neural cell adhesion molecule",
  "nmda receptor",
  "ampa receptor",
  "gaba receptor",
  "neural cell adhesion molecule l1",
  "microtubule-associated protein tau",
  "sphingomyelinase",
  "acid sphingomyelinase-like phosphodiesterase 3b"
)

rx_brain_ms <- regex(paste(brain_terms, collapse = "|"), ignore_case = TRUE)

# Create interest_MS flag
sig_ms <- sig_ms %>%
  mutate(
    interest_MS = dplyr::case_when(
      # EBV
      group_chr %in% ebv_ms_names | str_detect(group_lc, rx_ebv_ms) ~ "EBV",
      # CMV
      group_chr %in% cmv_ms_names | str_detect(group_lc, rx_cmv_ms) ~ "CMV",
      # Akkermansia
      group_chr %in% akk_ms_names | str_detect(group_lc, rx_akk_ms) ~ "Akkermansia",
      # Human brain-related antigens (if you have Origin_class == "human")
      (!("origin-class" %in% names(sig_ms)) |
         origin_class == "human") &
        (group_chr %in% human_brain_names | str_detect(group_lc, rx_brain_ms)) ~ "human_brain",
      TRUE ~ "none"
    )
  )

# Quick check
table(sig_ms$interest_MS)


#!/usr/bin/env Rscript

# Step 6: visualization of standardized residuals by antigen label
# Uses output from step 5: fisher_labeled_interestMS1.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridis)
})


# ================== ORIGINAL PLOTTING CODE (unchanged) ==================

# Prep
df_plot <- sig_ms %>%
  #filter(origin_class != "other", !is.na(std_residual)) %>%
  mutate(lev = factor(lev))

df_plot <- sig_ms %>%
  mutate(
    lev = factor(lev),
    origin_class = factor(
      origin_class,
      levels = c("human", "viral", "bacteria", "other")
    )
  )




df_hv   <- df_plot %>% filter(origin_class %in% c("human","viral","bacteria", "other"))
df_rest <- df_plot %>% filter(!origin_class %in% c("human","viral","bacteria", "other"))

# Plot: box trasparente (solo bordo) + jitter colorato per lev
ggplot(df_plot, aes(x = origin_class, y = std_residual)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, alpha = 0.6) +
  geom_boxplot(width = 0.65, outlier.shape = NA,
               fill = NA, color = "black", linewidth = 0.5) +  
   # jitter "normale" per le etichette meno dense
  geom_point(
    data = df_rest,
    aes(color = lev),
    position = position_jitter(width = 0.15, height = 0.08, seed = 123),
    size = 1.0,
    alpha = 0.75
  ) +
  # jitter più ampio (x e y) per human/viral
  geom_point(
    data = df_hv,
    aes(color = lev),
    position = position_jitter(width = 0.35, height = 0.15, seed = 456),
    size = 1.0,
    alpha = 0.65
  )  +
  scale_color_viridis_d(option = "D", end = 0.95) +
  labs(
    x = NULL,
    y = "Standardized residual (antigen enrichment)",
    color = "lev"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank()
  )

