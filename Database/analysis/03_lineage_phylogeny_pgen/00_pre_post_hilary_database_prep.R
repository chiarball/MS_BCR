#!/usr/bin/env Rscript

## 01_hilary_pre_post.R
## Pre-processing for HILARy (add sequence_id, filter bad rows)
## and post-processing of HILARy output to map clone_id ↔ cluster_h1 for selected clusters.

library(data.table)
library(dplyr)
library(stringr)

###############################################################################
## 1) PRE-HILARY: build MS_HC_naive_SEQID.tsv
###############################################################################

# Original naive DB (without sequence_id)
MS_HC_H <- fread("/doctorai/chiarba/analysis/MS_HC_naive.tsv")

unique(MS_HC_H$study)

# Remove unwanted studies
to_remove <- c("Stern", "Ramesh", "Palanichemy", "Greenfield", "Ghraichy")

n_before <- nrow(MS_HC_H)
MS_HC_H <- MS_HC_H[!tolower(trimws(study)) %in% tolower(to_remove)]
message("Removed ", n_before - nrow(MS_HC_H), " rows; remaining: ", nrow(MS_HC_H))

# Old file already prepared for HILARy (with sequence_id)
for_H_old <- fread("/doctorai/chiarba/bcr_from_raw/MS_HC_for_hilary.tsv")

# Ensure sequence_id is character
for_H_old <- for_H_old %>% mutate(sequence_id = as.character(sequence_id))
MS_HC_H   <- MS_HC_H   %>% mutate(sequence_id = as.character(sequence_id))

# Merge old + new
for_H_new <- bind_rows(for_H_old, MS_HC_H)

# Drop Ghraichy if present in the merged object
for_H_new <- for_H_new %>%
  filter(!(tolower(study) == "ghraichy"))

# Remove any existing seq_id / sequence_id columns to regenerate clean IDs
for_H_new[, (intersect(names(for_H_new), c("seq_id", "sequence_id"))) := NULL]

unique(for_H_new$study)

# Parameters for sequence_id creation
abbrev_len    <- 3   # number of letters from study to use as prefix
pad_digits    <- 7   # zero-padding digits (7 handles up to ~10M per study)
use_underscore <- TRUE

# Replace missing/empty study with "UNK"
for_H_new[is.na(study) | study == "", study := "UNK"]

# Abbreviation from study (remove non-alphanumeric, take first abbrev_len chars, upper case)
for_H_new[, .study_abbrev := toupper(substr(gsub("[^A-Za-z0-9]", "", study), 1, abbrev_len))]

# Per-study running number with zero padding
fmt <- paste0("%0", pad_digits, "d")
for_H_new[, .seq_num := sprintf(fmt, seq_len(.N)), by = .study_abbrev]

# Build sequence_id
if (use_underscore) {
  for_H_new[, sequence_id := paste0(.study_abbrev, "_", .seq_num)]
} else {
  for_H_new[, sequence_id := paste0(.study_abbrev, .seq_num)]
}

# Drop helper columns
for_H_new[, c(".study_abbrev", ".seq_num") := NULL]

print(head(for_H_new[, .(study, sequence_id)], 10))

# Check CDR3 / alignments and filter bad rows
sum(is.na(for_H_new$cdr3) | trimws(for_H_new$cdr3) == "")
for_H_new <- dplyr::filter(for_H_new, !is.na(cdr3) & trimws(cdr3) != "")

sum(is.na(for_H_new$sequence_alignment) | for_H_new$sequence_alignment == "")
sum(is.na(for_H_new$germline_alignment) | for_H_new$germline_alignment == "")

# Write final input for HILARy
fwrite(for_H_new,
       file = "/doctorai/chiarba/analysis/hilary/MS_HC_naive_SEQID.tsv",
       sep  = "\t")

## OPTIONAL: example for restricting to some subjects (keep commented)
## keep <- c("Ru8")
## for_H_last2 <- dplyr::filter(for_H_new, subject_number %in% keep)
## fwrite(for_H_last2,
##        file = "/doctorai/chiarba/analysis/hilary/MS_HC_sub_last2.tsv",
##        sep  = "\t")

###############################################################################
## 2) POST-HILARY: map clone_id ↔ cluster_h1 for selected clusters (RAxML input)
###############################################################################

# Paths
clusters_path <- "/doctorai/chiarba/analysis/MS_HC_naive_h0h1.tsv"
out_dir       <- "/doctorai/chiarba/analysis/hilary/output"
out_file      <- file.path(out_dir, "ALL_clone_ids_by_cluster_h1.tsv")

# List of cluster_h1 you want to use (unchanged)
selected_clusters <- c(
  520174,3944773,1332354,460740,1910226,2069857,3942486,2469355,
  579115,840865,3025717,361365,824794,3028615,3070060,3075043,
  3284237,484536,1926605,5889,27479,599754,2247404,3951452,
  448379,498087,3446945,767481,1364198,2414503,773274,2504107,
  2705643,11891,242486,251877,326178,577022,834074,1439725,
  1925874,1994679,3035626,205410,376210,519964,542137,557765,
  710219,817941,984205,1052183,1197602,1826283,2245214,2432527,
  2492934,2494325,3040943,3811523,34377,251876,274822,298656,
  334606,528889,570258,597204,648278,708689,764852,766091,
  778699,780366,814330,1018960,1600081,1753082,1754395,1764599,
  1789777,1842813,2386002,2565795,2607863,2661290,2676669,2706998,
  2716260,2808709,3021661,3146039,3184036,3342207,3701650,3739381,
  3810965,15556,221434,247819
)

# 2.1) Read clusters (with cluster_h1 + sequence)
clusters <- fread(clusters_path, sep = "\t", na.strings = c("", "NA"))

# Pick first available sequence-like column
seq_candidates <- c("sequence","v_cdr3_j","cdr3_aa","junction_aa","junction","cdr3")
seq_col_clusters <- seq_candidates[seq_candidates %in% names(clusters)][1]
if (is.na(seq_col_clusters)) stop("No sequence column found in clusters file.")

# Keep only selected clusters (comment this filter out if you want all)
clusters <- clusters[
  cluster_h1 %in% as.character(selected_clusters) |
  cluster_h1 %in% selected_clusters
]

# Standardize sequence for matching (trim + uppercase)
clusters_small <- unique(
  clusters[, .(
    sequence_std = toupper(trimws(get(seq_col_clusters))),
    cluster_h1   = as.character(cluster_h1)
  )]
)[!is.na(sequence_std) & sequence_std != ""]

# 2.2) Read and combine ALL HILARy outputs
files <- list.files(out_dir, pattern = "\\.tsv$", full.names = TRUE)
if (!length(files)) stop("No .tsv files found in output folder.")

hilary <- rbindlist(
  lapply(files, fread, sep = "\t", na.strings = c("", "NA")),
  fill = TRUE
)

# Detect sequence and clone_id columns in HILARy outputs
seq_col_hilary   <- seq_candidates[seq_candidates %in% names(hilary)][1]
clone_candidates <- c("clone_id","cloneid","clonotype_id","clonotype","clone")
clone_col_hilary <- clone_candidates[clone_candidates %in% names(hilary)][1]

if (is.na(seq_col_hilary))   stop("No sequence column found in HILARy outputs.")
if (is.na(clone_col_hilary)) stop("No clone_id-like column found in HILARy outputs.")

hilary_min <- hilary[, .(
  sequence_std = toupper(trimws(get(seq_col_hilary))),
  clone_id     = as.character(get(clone_col_hilary))
)][!is.na(sequence_std) & sequence_std != ""]

# 2.3) LEFT JOIN by sequence: clusters → hilary to get clone_id
merged <- merge(
  clusters_small, hilary_min,
  by = "sequence_std", all.x = TRUE, allow.cartesian = TRUE
)

# Keep only clone_id + cluster_h1, drop rows without clone_id
clones_db <- unique(merged[!is.na(clone_id), .(clone_id, cluster_h1)])

# 2.4) Save
fwrite(clones_db, out_file, sep = "\t", quote = FALSE, na = "NA")

# 2.5) Small summary
cat("Clusters used:   ", uniqueN(clusters_small$cluster_h1), "\n")
cat("Sequences used:  ", nrow(clusters_small), "\n")
cat("Clone IDs found: ", nrow(clones_db), "\n")
cat("Saved to:        ", out_file, "\n")
