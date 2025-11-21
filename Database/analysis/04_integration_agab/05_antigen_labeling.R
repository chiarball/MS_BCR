#!/usr/bin/env Rscript

# Step 5: classify antigen groups into broad origin classes (label)
# - Input:  fisher_results_annot1_antigen_all_levels_sig.csv
# - Output: fisher_labeled_interestMS1.tsv
# Code comments in English only.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
})

# ---------------------- I/O PATHS ----------------------

base_dir  <- "/doctorai/chiarba/AbAg_database/clean"
in_path   <- file.path(base_dir, "fisher_results_annot1_antigen_all_levels_sig.csv")
out_path  <- file.path(base_dir, "fisher_labeled_interestMS1.tsv")

# ---------------------- LOAD DATA ----------------------

sig_ms <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)

if (!"group" %in% names(sig_ms)) {
  stop("Column 'group' not found in input file. Check Fisher output.")
}

# ---------------------- BUILD GROUP TABLE ----------------------

# Unique group names to classify
gvec <- unique(as.character(sig_ms$group))
df   <- tibble(group = gvec, group_lc = tolower(gvec))

# ---------------------- REGEX RULES (MERGED VERSION) ----------------------

# 1) Toxins
rx_toxin <- c(
  "\\btoxin\\b", "neurotoxin", "exotoxin", "endotoxin",
  "alpha-?hemolysin", "aureolysin", "dermonecrotic", "vaginolysin", "botulinum"
)

# 2) Parasites (protozoa/helminths)
rx_parasite <- c(
  "plasmodium|merozoite|ookinete|circumsporozoite",
  "toxoplasma|\\bsag1\\b", "leishmania", "trypanosom|\\bvsg\\b",
  "schisto", "echinococcus", "taenia", "onchocerca", "giardia", "entamoeba", "trichomon"
)

# 3) Viral (more inclusive merged list)
rx_viral <- c(
  # Explicit Spike first
  "\\bspike glycoprotein\\b",
  # EBV / CMV
  "epstein[- ]?barr|\\bebv\\b",
  "cytomegalovirus|\\bhcmv\\b|\\bcmv\\b",
  # Other viruses
  "adenovirus|herpes(virus)?",
  "influenza|hemagglutinin|neuraminidase",
  "\\bhiv\\b|envelope glycoprotein gp\\d+|gp160",
  "papillomavirus|\\bhpv\\b",
  "hepatitis [abcdeg]|\\bhbv\\b|\\bhcv\\b",
  "coronavirus|sars[-_ ]?cov\\w*",
  # Generic viral ORFs and structural proteins
  "\\borf\\d+\\b|orf3a|orf7a|orf9b",
  "\\bnucleoprotein\\b|\\bcapsid protein\\b|\\benvelope protein\\b|\\bmatrix protein\\b",
  "\\breplicase\\b|polyprotein",
  # AAV family
  "\\baav\\b|adeno[- ]associated virus",
  "assembly activating protein|\\baap\\b",
  "protein rep\\s?(40|52|68|78)\\b|\\brep(40|52|68|78)\\b|\\bvp1\\b"
)

# 4) Bacterial (merged list, includes Akkermansia)
rx_bacterial <- c(
  "akkermansia", "staphyl", "streptoc", "pseudomonas", "mycobacter", "clostrid",
  "escherichia|\\becoli\\b", "salmonella", "listeria", "neisseria", "borrelia",
  "brucella", "vibrio", "campylobacter", "bacteroides",
  "outer membrane", "lipopolysaccharide|\\blps\\b"
)

# 5) Plant
rx_plant <- c(
  "pollen|allergen", "\\bphl p\\b", "bet v 1", "\\blol p\\b",
  "arabidopsis", "\\bplant\\b"
)

# 6) Human (core) – kept broad but with viral exclusions
rx_human_core <- c(
  "\\bhla\\b|histocompatibility",
  "\\bcd\\d+\\b",
  "interleukin|\\bil-?\\d",
  "tumor necrosis factor|\\btnf\\b",
  "chemokine|\\b(?:cxcl|ccl|xcl|cx3cl)\\d+\\b",
  "immunoglobulin|heavy constant|kappa\\b|lambda\\b|fc receptor|fcrn|mhc",
  "receptor|integrin|kinase|phosphatase|transporter|enzyme|growth factor|hormone|collagen|matrix|protease|glycoprotein",
  "neurotensin receptor|adrenergic|dopamine|serotonin|acetylcholine receptor",
  "albumin|ferritin|transferrin|actinin|myosin|dystroph|sarcoglycan|utrophin|laminin",
  "beta-?amyloid|\\btau\\b|alpha-?synuclein|prion|psen\\d|pink1",
  "insulin|leptin|ghrelin|calcitonin|glucagon|thyroid|prolactin",
  "\\bjak\\b|\\bstat\\b|mapk|erk|\\bakt\\b|raf\\b|pi3k|mtor|tgf-?beta|bmp\\b|smad\\b|\\bwnt\\b|frizzled|\\bnotch\\b",
  "carboxylesterase 1\\b|\\blysozyme c\\b"
)

# Words that, if present, should block a "human" label
viral_exclusion <- "(spike glycoprotein|capsid|nucleoprotein|\\borf\\d+\\b|polyprotein|replicase|aav|gp\\d+)"

# ---------------------- CLASSIFIER ----------------------

classify_group <- function(x) {
  # x is lowercased group name
  pick_any <- function(rx_vec) {
    any(str_detect(x, regex(paste(rx_vec, collapse = "|"), ignore_case = TRUE)))
  }
  if (pick_any(rx_toxin))     return("toxin")
  if (pick_any(rx_parasite))  return("parasite")
  if (pick_any(rx_viral))     return("viral")
  if (pick_any(rx_bacterial)) return("bacterial")
  # human only if no clear viral markers
  if (!str_detect(x, regex(viral_exclusion, ignore_case = TRUE)) &&
      pick_any(rx_human_core)) return("human")
  if (pick_any(rx_plant))     return("plant")
  return("other")
}

df <- df %>%
  mutate(broad_origin = vapply(group_lc, classify_group, character(1)))

# ---------------------- TARGETED OVERRIDES ----------------------

override_map <- c(
  # explicit fixes
  "lysozyme c"         = "human",
  "spike glycoprotein" = "viral",
  "protein rep52"      = "viral",
  # EBV/CMV/Akkermansia (reinforce classification)
  "epstein-barr"       = "viral",
  "ebv"                = "viral",
  "cytomegalovirus"    = "viral",
  "akkermansia"        = "bacterial",
  # human enzymes that may look generic
  "liver carboxylesterase 1" = "human"
)

for (k in names(override_map)) {
  df$broad_origin[str_detect(df$group_lc, fixed(k, ignore_case = TRUE))] <- override_map[[k]]
}

# ---------------------- JOIN BACK TO sig_ms ----------------------

sig_ms <- sig_ms %>%
  mutate(group_chr = as.character(group)) %>%
  left_join(select(df, group, broad_origin), by = c("group_chr" = "group")) %>%
  select(-group_chr)

# For downstream plotting we will use 'label' (same as broad_origin)
sig_ms <- sig_ms %>%
  mutate(label = broad_origin)

# ---------------------- QUICK SUMMARY ----------------------

cat("Unique groups per class:\n")
print(table(df$broad_origin))

cat("\nCounts on full sig_ms rows:\n")
print(
  sig_ms %>%
    count(label, name = "n_sequences") %>%
    arrange(desc(n_sequences))
)

# ---------------------- SAVE OUTPUT ----------------------

write_tsv(sig_ms, out_path)
cat("\nSaved labeled table to:", out_path, "\n")

