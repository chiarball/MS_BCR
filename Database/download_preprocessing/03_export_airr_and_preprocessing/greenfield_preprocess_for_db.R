#!/usr/bin/env Rscript

# greenfield_preprocess_for_db.R
# Preprocessing of Greenfield (baseline and T2) for DB.

library(data.table)
library(dplyr)
library(stringr)
library(tools)

#### 1) BASELINE (Green_all.tsv) #############################################

db <- fread("/doctorai/chiarba/bcr_from_raw/green/AIRR/Green_all.tsv")

metadata <- fread("/doctorai/chiarba/bcr_from_raw/green/AIRR/Metadata_green_SraRunTable.csv")

metadata <- metadata %>%
  rename(SRA_id = Run)

# Select useful columns from metadata
meta_subset <- metadata %>%
  select(
    SRA_id,
    sex,
    tissue,
    age,
    `Library Name`,
    `Sample Name`,
    BioProject,
    BioSample,
    cell_subtype,
    disease
  ) %>%
  rename(AGE = age)

# Merge AIRR table with metadata
db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")

# Remove TCR annotations (keep only BCR)
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

# subject_number from "Library Name": PtX -> GrX
db_meta <- db_meta %>%
  rename(subject_number = `Library Name`) %>%
  mutate(subject_number = str_extract(subject_number, "^Pt\\d+")) %>%
  mutate(subject_number = str_replace(subject_number, "^Pt", "Gr"))

# study label
db_meta <- db_meta %>%
  mutate(study = "Greenfield")

# Harmonize tissue
db_meta <- db_meta %>%
  mutate(tissue = case_when(
    tissue == "cerebrospinal fluid" ~ "CSF",
    tissue == "peripheral blood"    ~ "blood",
    TRUE                            ~ tissue
  ))

# disease_diagnosis
db_meta <- db_meta %>%
  rename(disease_diagnosis = disease) %>%
  mutate(disease_diagnosis = case_when(
    disease_diagnosis == "multiple sclerosis" ~ "MS",
    TRUE                                      ~ disease_diagnosis
  ))

# Preserve original sequence_id
db_meta <- rename(db_meta, sequence_id_old = sequence_id)

# Write baseline DB-ready table
fwrite(
  db_meta,
  file = "/doctorai/chiarba/dataset/green_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)


#### 2) T2 (green_add_T2/AIRR) ###############################################

folder_path_t2 <- "/doctorai/chiarba/bcr_from_raw/green/green_add_T2/AIRR"
tsv_files_t2   <- list.files(folder_path_t2, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- file_path_sans_ext(basename(file_path))
  df <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

# Same logic as map_dfr: read each file and row-bind into a single data.table
db_t2_list <- lapply(tsv_files_t2, read_and_tag)
db_t2 <- rbindlist(db_t2_list, use.names = TRUE, fill = TRUE)

metadata_t2 <- fread("/doctorai/chiarba/bcr_from_raw/green/green_add_T2/metadata_greenT2.csv") %>%
  rename(SRA_id = Run)

# Select useful columns from metadata
meta_subset_t2 <- metadata_t2 %>%
  select(
    SRA_id,
    sex,
    tissue,
    AGE,
    `Library Name`,
    `Sample Name`,
    BioProject,
    BioSample,
    cell_subtype,
    disease
  )

# Merge AIRR table with metadata
db_meta_t2 <- db_t2 %>%
  left_join(meta_subset_t2, by = "SRA_id")

# Remove TCR annotations (keep only BCR)
any(grepl("^TR", db_meta_t2$v_call))
any(grepl("^TR", db_meta_t2$d_call))
any(grepl("^TR", db_meta_t2$j_call))

db_meta_t2 <- db_meta_t2[
  !(
    grepl("^TR", ifelse(is.na(v_call), "", v_call), ignore.case = TRUE) |
    grepl("^TR", ifelse(is.na(d_call), "", d_call), ignore.case = TRUE) |
    grepl("^TR", ifelse(is.na(j_call), "", j_call), ignore.case = TRUE)
  )
]

# subject_number from "Library Name": PtX -> Gr_T2_X
db_meta_t2 <- db_meta_t2 %>%
  rename(subject_number = `Library Name`) %>%
  mutate(subject_number = str_extract(subject_number, "^Pt\\d+")) %>%
  mutate(subject_number = str_replace(subject_number, "^Pt", "Gr_T2_"))

# study label
db_meta_t2 <- db_meta_t2 %>%
  mutate(study = "Greenfield")

# Harmonize tissue
db_meta_t2 <- db_meta_t2 %>%
  mutate(tissue = case_when(
    tissue == "cerebrospinal fluid" ~ "CSF",
    tissue == "peripheral blood"    ~ "blood",
    TRUE                            ~ tissue
  ))

# disease_diagnosis
db_meta_t2 <- db_meta_t2 %>%
  rename(disease_diagnosis = disease) %>%
  mutate(disease_diagnosis = case_when(
    disease_diagnosis == "multiple sclerosis" ~ "MS",
    TRUE                                      ~ disease_diagnosis
  ))

# Preserve original sequence_id
db_meta_t2 <- rename(db_meta_t2, sequence_id_old = sequence_id)

# Write T2 DB-ready table
fwrite(
  db_meta_t2,
  file = "/doctorai/chiarba/dataset/green_for_db2.tsv",
  sep = "\t",
  row.names = FALSE
)
