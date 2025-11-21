#!/usr/bin/env Rscript

# CDR3 extraction from AgAb database (clean version)

library(data.table)
library(dplyr)
library(stringr)

# Input and output paths
base_dir <- "/doctorai/chiarba/AbAg_database"
out_dir  <- file.path(base_dir, "clean")
if (!dir.exists(out_dir)) dir.create(out_dir)

agab <- fread(file.path(base_dir, "agab_merged.csv"))

# Allowed AA and motif patterns
aa_stdX   <- "ACDEFGHIKLMNPQRSTVWYX"
const_pat <- "(?:VTVS[ST]|SAS?TK)"
j_pat     <- paste0("W[", aa_stdX, "]{1,3}G")
start_pat <- "(?:YYC|YFC|FYC|NYC|HYC|RFC)"

# Extract CDR3 from a single sequence
extract_cdr3_one_v5 <- function(s) {
  if (is.na(s) || nchar(s) < 20) return(NA_character_)
  s <- toupper(gsub("[^A-Z]", "", s))

  j_all <- gregexpr(j_pat, s, perl = TRUE)[[1]]
  if (j_all[1] == -1) return(NA_character_)

  const_m <- regexpr(const_pat, s, perl = TRUE)[1]
  if (const_m != -1) {
    j_all <- as.integer(j_all)
    j_before_const <- j_all[j_all < const_m]
    if (length(j_before_const) == 0) return(NA_character_)
    j_start <- max(j_before_const)
  } else {
    j_start <- max(as.integer(j_all))
  }

  pre_j <- if (j_start > 1L) substr(s, 1, j_start - 1L) else ""
  a_hits <- gregexpr(start_pat, pre_j, perl = TRUE)[[1]]

  if (a_hits[1] != -1) {
    c_pos <- max(as.integer(a_hits)) + 2L
  } else {
    c_hits <- gregexpr("C", pre_j, perl = TRUE)[[1]]
    if (c_hits[1] == -1) return(NA_character_)
    c_positions <- as.integer(c_hits)
    lens <- (j_start - 1L) - c_positions + 1L
    ok <- lens >= 3L & lens <= 80L
    if (!any(ok)) return(NA_character_)
    c_pos <- max(c_positions[ok])
  }

  end_pos <- j_start - 1L
  if (end_pos < c_pos) return(NA_character_)

  cdr3 <- substr(s, c_pos, end_pos)
  if (!grepl(paste0("^[", aa_stdX, "]+$"), cdr3)) return(NA_character_)
  cdr3
}

# Vectorized wrapper
extract_cdr3_vec_v5 <- function(v) vapply(v, extract_cdr3_one_v5, character(1))

# Apply extraction
agab <- agab %>%
  mutate(
    heavy_sequence = toupper(heavy_sequence),
    cdr3 = extract_cdr3_vec_v5(heavy_sequence),
    cdr3_len = nchar(cdr3),
    cdr3_filtered = if_else(!is.na(cdr3) & cdr3_len >= 5 & cdr3_len <= 30, cdr3, NA_character_)
  )

# Remove leading C and trailing W, rename output columns
agab <- agab %>%
  mutate(
    cdr3 = if_else(!is.na(cdr3),
                   str_replace(str_replace(cdr3, "^C", ""), "W$", ""),
                   NA_character_),
    cdr3_len = nchar(cdr3)
  ) %>%
  rename(
    cdr3_aa = cdr3,
    cdr3_aa_len = cdr3_len
  )

# Save output
fwrite(agab, file.path(out_dir, "agab_cdr3_extracted1.tsv"), sep = "\t", row.names = FALSE)
