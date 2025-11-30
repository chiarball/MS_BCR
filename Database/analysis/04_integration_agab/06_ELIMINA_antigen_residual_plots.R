#!/usr/bin/env Rscript

# Step 6: visualization of standardized residuals by antigen label
# Uses output from step 5: fisher_labeled_interestMS1.tsv
# Comments in English only.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridis)
})

# ---- I/O ----
base_dir <- "/doctorai/chiarba/AbAg_database/clean"
in_path  <- file.path(base_dir, "fisher_labeled_interestMS1.tsv")

sig_labeled <- read_tsv(in_path, show_col_types = FALSE, guess_max = 100000)

# Prep
df_plot <- sig_labeled %>%
  filter(label != "other", !is.na(std_residual)) %>%
  mutate(lev = factor(lev))

# Order x by decreasing count
ord_x <- df_plot %>%
  count(label, name = "n") %>%
  arrange(desc(n)) %>%
  pull(label)

df_plot <- df_plot %>%
  mutate(label = factor(label, levels = ord_x))

df_hv   <- df_plot %>% filter(label %in% c("human","viral"))
df_rest <- df_plot %>% filter(!label %in% c("human","viral"))

ggplot(df_plot, aes(x = label, y = std_residual)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, alpha = 0.6) +
  geom_boxplot(width = 0.65, outlier.shape = NA,
               fill = NA, color = "black", linewidth = 0.5) +
  geom_point(
    data = df_rest,
    aes(color = lev),
    position = position_jitter(width = 0.15, height = 0.08, seed = 123),
    size = 1.1,
    alpha = 0.75
  ) +
  geom_point(
    data = df_hv,
    aes(color = lev),
    position = position_jitter(width = 0.35, height = 0.15, seed = 456),
    size = 1.1,
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
