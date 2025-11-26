#!/usr/bin/env Rscript

# agrafiotis_preprocess_for_db.R
# Preprocessing of Agrafiotis single-cell data.

library(data.table)
library(dplyr)
library(stringr)

## 1. Merge AIRR-format contig files and metadata 

path_airr  <- "/doctorai/chiarba/add_db/agrafiotis/AIRR_format"
files_airr <- list.files(path_airr, pattern = "\\.tsv$", full.names = TRUE)

agrafiotis <- rbindlist(
  lapply(files_airr, function(f) {
    dt <- fread(f, sep = "\t", na.strings = c("", "NA"))
    dt[, `__source_file` := basename(f)]
    dt
  }),
  fill = TRUE
)

metadata_arg <- fread("/doctorai/chiarba/add_db/agrafiotis/metadataagrafiotis.csv")

metadata_arg <- metadata_arg %>%
  rename(
    SRA_id   = Run,
    sample_id = `Sample Name`
  )

# (Facoltativo, non usato dopo ma non fa danni)
agrafiotis$SRR <- sub(
  ".*_(SRR[0-9]+)\\.contigs\\.tsv$",
  "\\1",
  agrafiotis$`__source_file`
)

# Extract SRA_id from file name (SRRxxxxx)
agrafiotis <- agrafiotis %>%
  mutate(
    SRA_id = str_extract(`__source_file`, "SRR\\d+")
  )

# Select useful metadata columns
meta_subset1 <- metadata_arg %>%
  select(
    SRA_id,
    sex,
    tissue,
    AGE,
    disease_stage,
    sample_id,
    disease,
    BioSample,
    BioProject
  )

agr_sra_meta <- agrafiotis %>%
  left_join(meta_subset1, by = "SRA_id")

## 2. Single-cell cleaning: filter by cell_id and collapse by clone_id 

# Keep only cell_id with 1 or 2 rows per SRA_id–cell_id
agr_sra_meta <- agr_sra_meta %>%
  group_by(SRA_id, cell_id) %>%
  mutate(n_per_cell = n()) %>%
  ungroup() %>%
  filter(n_per_cell <= 2) %>%
  select(-n_per_cell)

# Collapse to one row per SRA_id–clone_id (keep first)
agr_collapsed <- agr_sra_meta %>%
  group_by(SRA_id, clone_id) %>%
  mutate(n_cells = n()) %>%
  ungroup() %>%
  distinct(SRA_id, clone_id, .keep_all = TRUE)

fwrite(
  agr_collapsed,
  "/doctorai/chiarba/add_db/agrafiotis/agr_cellid_cloneid.tsv",
  sep = "\t",
  row.names = FALSE
)

## 3. Preprocessing for DB 

# Check for TCR annotations
any(grepl("^TR", agr_collapsed$v_call))
any(grepl("^TR", agr_collapsed$d_call))
any(grepl("^TR", agr_collapsed$j_call))

# Remove TR in V/D/J
agr_collapsed <- agr_collapsed %>%
  filter(
    !str_detect(ifelse(is.na(d_call), "", d_call), "^TR"),
    !str_detect(ifelse(is.na(v_call), "", v_call), "^TR"),
    !str_detect(ifelse(is.na(j_call), "", j_call), "^TR")
  )

# subject_number: "X-123" -> "ArX123" using sample_id
agr_collapsed <- agr_collapsed %>%
  mutate(
    subject_number = paste0(
      "Ar",
      substr(sample_id, 1, 1),
      sub(".*-(\\d+)$", "\\1", sample_id)
    )
  )

# study column
agr_collapsed <- agr_collapsed %>%
  mutate(study = "Agrafiotis")

# disease_diagnosis from disease_stage 
agr_collapsed <- agr_collapsed %>%
  rename(disease_diagnosis = disease_stage)

# Preserve original sequence_id
agr_collapsed <- rename(agr_collapsed, sequence_id_old = sequence_id)

# Write final DB-ready Agrafiotis 
fwrite(
  agr_collapsed,
  "/doctorai/chiarba/dataset/agr_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)

