#!/usr/bin/env Rscript

# greenfield_preprocess_for_db.R
# Preprocessing of Greenfield for DB.

library(data.table)
library(dplyr)
library(stringr)

# 1. Load AIRR data and metadata

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
  )

meta_subset <- meta_subset %>%
  rename(AGE = age)

# Merge AIRR table with metadata
db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")

# 2. Remove TCR annotations (keep only BCR)

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

# 3. Create subject_number and study columns

# subject_number from "Library Name": PtX -> GrX
db_meta <- db_meta %>%
  rename(subject_number = `Library Name`) %>%
  mutate(subject_number = str_extract(subject_number, "^Pt\\d+"))

db_meta <- db_meta %>%
  mutate(subject_number = str_replace(subject_number, "^Pt", "Gr"))

# study label
db_meta <- db_meta %>%
  mutate(study = "Greenfield")

# 4. Harmonize tissue and disease_diagnosis annotations

# tissue
db_meta <- db_meta %>%
  mutate(tissue = case_when(
    tissue == "cerebrospinal fluid" ~ "CSF",
    tissue == "peripheral blood"    ~ "blood",
    TRUE                            ~ tissue   # keep other values unchanged
  ))

# disease_diagnosis
db_meta <- db_meta %>%
  rename(disease_diagnosis = disease) %>%
  mutate(disease_diagnosis = case_when(
    disease_diagnosis == "multiple sclerosis" ~ "MS",
    TRUE                                      ~ disease_diagnosis
  ))

# 5. Preserve original sequence_id and write out DB-ready table

db_meta <- rename(db_meta, sequence_id_old = sequence_id)

fwrite(
  db_meta,
  file = "/doctorai/chiarba/dataset/green_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)
