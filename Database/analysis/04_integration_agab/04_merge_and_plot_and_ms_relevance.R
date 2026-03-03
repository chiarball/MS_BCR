#!/usr/bin/env Rscript

# Merge Fisher antigen-level results (lev0/1/2) and keep only significant rows
# Comments in English only

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# --- INPUT BASE DIRECTORY (clean version) ---
base_dir <- "/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test"

# --- INPUT FILES (from step 3 / Python Fisher) ---
files <- c(
  lev0 = file.path(base_dir, "fisher_results_annot2_antigen_lev0.csv"),
  lev1 = file.path(base_dir, "fisher_results_annot2_antigen_lev1.csv"),
  lev2 = file.path(base_dir, "fisher_results_annot2_antigen_lev2.csv")
)

# --- READ EACH LEVEL AND ADD 'lev' COLUMN ---
dfs <- lapply(names(files), function(nm) {
  df <- fread(files[[nm]])
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
fwrite(sig, out_path)

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

sig_ms <- fread("/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test/fisher_results_annot2_antigen_all_levels_sig.csv")

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
  library(ggplot2, lib.loc ="/doctorai/niccoloc/libR2")
  library(labeling, lib.loc ="/doctorai/niccoloc/libR2")
  library(farver, lib.loc ="/doctorai/niccoloc/libR2")
  library(ggrepel, lib.loc ="/doctorai/niccoloc/libR2")
  library(viridis)
})

#manual annotation for crashed entries
origin_manual <- tibble::tribble(
  ~group_lc, ~origin_class,
  "lysozyme c", "human",
  "spike glycoprotein", "viral",
  "immunoglobulin kappa constant", "human",
  "superfamily_necrosis_factor_member_tumor_hepatitis_cellular_ligand", "human",
  "tissue factor", "human",
  "interleukin_subunit_beta_integrin_alpha", "human",
  "interleukin_subunit_alpha_granulocyte_stimulating_macrophage_colony_factor", "human",
  "lymphocyte_antigen_binding_lectin_sialic_acid_cd19_like", "human",
  "interleukin_necrosis_factor_tumor_17a", "human",
  "ubiquitin carboxyl-terminal hydrolase 30", "human",
  "beta-klotho", "human",
  "immunoglobulin iota chain", "human",
  "angiopoietin_endothelial_vascular_factor_growth_form_long", "human",
  "genome polyprotein", "viral",
  "cytotoxic t-lymphocyte protein 4", "human",
  "mesothelin", "human",
  "hemagglutinin", "viral",
  "inactive tyrosine-protein kinase transmembrane receptor ror1", "human",
  "integrin_cathepsin_alpha_beta", "human",
  "epstein-barr nuclear antigen 1", "viral",
  "cd27 antigen", "human",
  "thrombopoietin", "human",
  "c-x-c chemokine receptor type 2", "human",
  "major prion protein", "human",
  "neurogenic locus notch homolog protein 1", "human",
  "coagulation_factor_prothrombin_viii", "human",
  "interferon gamma", "human",
  "programmed cell death 1 ligand 1", "human",
  "transmembrane protease serine 2", "human"
)


### Plotting

# fwrite(df_plot, "/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test/old_df_plot_result.csv")  
# df_old= fread("/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test/old_df_plot_result.csv")

df_plot <- df_old

df_plot <- sig_ms %>%
  mutate(
    lev = factor(lev),  # make sure lev is a factor
    origin_class= str_remove(origin_class, "\\|other"),
    origin_class= ifelse(
      str_detect(origin_class, "\\|viral"), "viral", origin_class
    ),
    # origin_class = factor(
    #   origin_class,
    #   levels = c("human", "viral", "bacteria", "other")
    # )
  )  %>%
  # ensure group_lc exists and is clean
  mutate(group_lc = tolower(trimws(group))) %>%
  left_join(origin_manual, by = "group_lc") %>% 
  mutate(
    origin_class = ifelse(
      is.na(origin_class.x),
      origin_class.y,
      origin_class.x
    )
  ) %>% select(-origin_class.x, -origin_class.y) %>% 
  mutate(
    lev = factor(
      lev,
      levels = c("lev0", "lev1", "lev2"),
      labels = c(
        "0",
        "1",
        "2"
      ) ) ) 


df_plot %>% filter(is.na(origin_class)) %>% pull(group) %>% unique()

df_hv   <- df_plot %>% filter(origin_class %in% c("human"))
df_rest <- df_plot %>% filter(origin_class %in% c("viral","bacteria", "other"))

p <- ggplot(df_plot  , aes(x = origin_class, y = std_residual)) +
  geom_boxplot(
    width = 0.65, outlier.shape = NA,
    fill = NA, color = "black", linewidth = 0.5
  ) +
  
  
  geom_point(
    data = df_plot %>% filter(interest_MS == "none"),
    aes(color = lev),
    position = position_jitter(width = 0.22, height = 0.12, seed = 123),
    size = 1.2,
    # star for ebv the rest normal
    
    alpha = 0.70,
    # show.legend = FALSE
  ) +
  # star for ebv the rest normal
  geom_point(
    data = df_plot %>% filter(interest_MS == "EBV"),
    aes(  shape = interest_MS),
    color= "52d8da",
 
    position = position_jitter(width = 0.22, height = 0.12, seed = 123),
    size = 5,
    # star for ebv the rest normal
    
    alpha = 1.2 
  ) +
  scale_shape_manual(
    name="",
    values = c(
      "EBV" = "★"        # This is the 'star' (asterisk)
 
     ))+
 

  scale_x_discrete(labels = c(
    human    = "HUMAN",
    viral    = "VIRUS",
    bacteria = "BACTERIA",
    other    = "OTHER"
  )) +
  
  labs(
    x = NULL,
    y = "Standardized residual (antigen enrichment)",
    color = "Levensthein\ndistance"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_blank()
  ) +
  
  scale_y_log10() +
  
  # facet_wrap(
  #   ~ lev,
  #   nrow = 1,
  #   labeller = labeller(
  #     lev = c(
  #       "lev0" = "Levenshtein distance = 0",
  #       "lev1" = "Levenshtein distance = 1",
  #       "lev2" = "Levenshtein distance = 2"
  #     )
  #   )
  # ) +
  
  # # mediane
  # stat_summary(
  #   fun = median,
  #   geom = "text",
  #   aes(label = sprintf("%.2f", after_stat(y))),
  #   size = 3.0,
  #   # color = "white",
  #   fontface = "bold",
  #   vjust = -0.6
  # ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = sprintf("%.2f", after_stat(y))),
    size = 2.0,
    color = "black",
    fontface = "bold",
    vjust = -0.6
  )  

p

grDevices::svg(
  filename = "/doctorai/niccoloc/agab.svg",
  width    = 9,
  height   = 7
)
print(p)

dev.off()


table(df_plot$origin_class)
