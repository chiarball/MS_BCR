#!/usr/bin/env Rscript

# laurent_preprocess_for_db.R
# Preprocessing of Laurent data for DB.

library(data.table)
library(dplyr)
library(stringr)

# 1. Read and merge all AIRR TSV files, tagging with SRA_id

folder_path <- "/doctorai/chiarba/add_db/Laurent/AIRR"
tsv_files   <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

db <- rbindlist(lapply(tsv_files, read_and_tag), use.names = TRUE, fill = TRUE)

# 2. Read metadata and select useful columns

metadata <- fread("/doctorai/chiarba/add_db/Laurent/metadata_lau.csv")

metadata <- metadata %>%
  rename(SRA_id = Run)

meta_subset <- metadata %>%
  select(
    SRA_id,
    sex,
    tissue,
    AGE,
    `Sample Name`,
    Patient_ID,
    receptor_type,
    time_point,
    treatment,
    weeks_since_baseline,
    BioSample,
    BioProject
  )

# 3. Merge AIRR table with metadata

db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")

# 4. Remove TCR annotations (keep only BCR)

any(grepl("^TR", db_meta$v_call))
any(grepl("^TR", db_meta$d_call))
any(grepl("^TR", db_meta$j_call))

db_meta <- db_meta[
  !(
    grepl("^TR", ifelse(is.na(v_call), "", v_call), ignore.case = TRUE) |
    grepl("^TR", ifelse(is.na(d_call), "", d_call), ignore.case = TRUE) |
    grepl("^TR", ifelse(is.na(j_call), "", j_call), ignore.case = TRUE)
  )
]

# 5. Create subject_number and study columns

db_meta <- db_meta %>%
  mutate(
    subject_number = paste0(
      "Lau",
      sub("MS-(\\d+)", "\\1", Patient_ID)
    )
  )

db_meta <- db_meta %>%
  mutate(study = "Laurent")

# 6. Harmonize tissue and disease_diagnosis

db_meta <- db_meta %>%
  mutate(tissue = "blood")

db_meta <- db_meta %>%
  mutate(disease_diagnosis = "MS")

# 7. Write out DB-ready laurent

db_meta <- rename(db_meta, sequence_id_old = sequence_id)

fwrite(
  db_meta,
  "/doctorai/chiarba/dataset/lau_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)
