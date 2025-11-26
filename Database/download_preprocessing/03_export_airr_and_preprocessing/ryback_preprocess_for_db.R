#!/usr/bin/env Rscript

# ryback_preprocess_for_db.R
# Preprocessing of Ryback AIRR dataset for unified DB.

library(data.table)
library(dplyr)
library(stringr)
library(tools)

# 1. Load and merge AIRR files 

folder_path <- "/doctorai/chiarba/add_db/Ryback/AIRR"

tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- file_path_sans_ext(basename(file_path))
  df <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

db <- rbindlist(lapply(tsv_files, read_and_tag), use.names = TRUE, fill = TRUE)

# 2. Load metadata and merge 

metadata <- fread("/doctorai/chiarba/add_db/Ryback/metadata_ryback.csv")
sex_age  <- fread("/doctorai/chiarba/add_db/Ryback/cohort_characteristics_SBL.csv")

setnames(metadata, "Submitter_Id", "sbl_id")
setnames(sex_age,   "SBL ID",      "sbl_id")

metadata <- merge(
  metadata,
  sex_age[, .(sbl_id, sex, agegroup)],
  by = "sbl_id",
  all.x = TRUE
)

metadata <- metadata %>%
  rename(SRA_id = Run)

meta_subset <- metadata %>%
  select(
    SRA_id,
    tissue_type,
    `Sample Name`,
    `Library Name`,
    sex,
    agegroup,
    BioProject,
    BioSample
  )

db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")


# 3. Remove TCR annotations 
# Some Ryback files use d_gene instead of d_call → handle both safely.

db_meta <- db_meta %>%
  filter(
    !str_detect(ifelse(is.na(v_call), "", v_call), "^TR"),
    !str_detect(ifelse(is.na(j_call), "", j_call), "^TR"),
    !str_detect(
      ifelse(is.na(ifelse("d_call" %in% names(db_meta), d_call, d_gene)), "",
             ifelse("d_call" %in% names(db_meta), d_call, d_gene)),
      "^TR"
    )
  )


# 4. Construct subject_number, harmonize fields 

db_meta <- db_meta %>%
  mutate(
    subject_number = if_else(
      !is.na(`Library Name`) & str_detect(`Library Name`, "SBL\\d+"),
      paste0("Ry", as.integer(str_extract(`Library Name`, "(?<=SBL)\\d+"))),
      NA_character_
    )
  )

db_meta <- db_meta %>%
  mutate(
    disease_diagnosis = "MS",
    tissue            = "blood",
    study             = "Ryback"
  )


# 5. Preserve original sequence_id and write out

db_meta <- rename(db_meta, sequence_id_old = sequence_id)

fwrite(
  db_meta,
  "/doctorai/chiarba/dataset/ryback_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)


