#!/usr/bin/env Rscript

# ruschill_preprocess_for_db.R
# Preprocessing of Ruschill for unified DB.

library(data.table)
library(dplyr)
library(stringr)
library(tools)

# 1. Load and merge AIRR files, tagging with SRA_id 

folder_path <- "/doctorai/chiarba/add_db/Ruschil/AIRR"

tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- file_path_sans_ext(basename(file_path))
  df <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

db <- rbindlist(lapply(tsv_files, read_and_tag), use.names = TRUE, fill = TRUE)

# 2. Load metadata and merge 

metadata <- fread("/doctorai/chiarba/add_db/Ruschil/Ruschil_metadata.csv")

metadata <- metadata %>%
  rename(SRA_id = Run)

meta_subset <- metadata %>%
  select(
    SRA_id,
    sex,
    tissue,
    AGE,
    `Sample Name`,
    BioProject,
    BioSample,
    cell_subtype,
    disease
  )

db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")

# 3. (Optional check) look for TCR annotations 

any(grepl("^TR", db_meta$v_call))
any(grepl("^TR", db_meta$d_call))
any(grepl("^TR", db_meta$j_call))
# N.B.: as in your original script, we only check, we do not filter TR here.

# 4. Add study, subject_number, harmonize tissue & diagnosis 

# study column
db_meta <- db_meta %>%
  mutate(study = "Ruschill")

# subject_number: from "Sample Name" of the form "...CLADXX..."
db_meta <- db_meta %>%
  mutate(
    subject_number = str_extract(`Sample Name`, "(?<=CLAD)\\d+") %>%
      paste0("Ru", .)
  )

# tissue harmonization
db_meta <- db_meta %>%
  mutate(tissue = case_when(
    tissue == "cerebrospinal fluid" ~ "CSF",
    tissue == "PBMC"               ~ "blood",
    TRUE                           ~ tissue
  ))

# disease_diagnosis
db_meta <- db_meta %>%
  mutate(disease_diagnosis = "MS")

# 5. Preserve original sequence_id and write out 

db_meta <- rename(db_meta, sequence_id_old = sequence_id)

fwrite(
  db_meta,
  "/doctorai/chiarba/dataset/rus_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)
