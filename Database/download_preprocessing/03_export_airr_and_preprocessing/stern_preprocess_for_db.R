#!/usr/bin/env Rscript

# stern_preprocess_for_db.R
# Preprocessing of Stern for DB.

library(data.table)
library(dplyr)
library(stringr)

# 1. Load AIRR data and metadata

db <- fread("/doctorai/chiarba/bcr_from_raw/stern/fastq_download/AIRR/stern_sra_db.tsv")

metadata <- fread("/doctorai/chiarba/bcr_from_raw/stern/meta_stern.csv")

metadata <- metadata %>%
  rename(SRA_id = Run)

# Select useful metadata
meta_subset <- metadata %>%
  select(
    SRA_id,
    `Library Name`,
    BioProject,
    BioSample,
    health_state
  )

# Merge AIRR table with metadata
db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")

# 2. Remove any TCR annotations (keep only BCR)

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

db_meta <- db_meta %>%
  mutate(
    subject_number = if_else(
      !is.na(`Sample Name`) & str_detect(`Sample Name`, "^M\\d+"),
      paste0("St", str_extract(`Sample Name`, "(?<=^M)\\d+")),
      NA_character_
    )
  )

db_meta <- db_meta %>%
  mutate(study = "Stern")

# 4. Harmonize tissue annotation

db_meta <- db_meta %>%
  mutate(tissue = case_when(
    tissue == "Cervical lymph node" ~ "CLN",
    TRUE                            ~ tissue
  ))

# 5. Create disease_diagnosis column from disease_stage

db_meta <- db_meta %>%
  mutate(
    disease_diagnosis = case_when(
      disease_stage == "Primary progressive"                  ~ "PPMS",
      disease_stage %in% c("Chronic progressive", "Chronic")  ~ "Chronic",
      disease_stage == "Secondary progressive"                ~ "SPMS",
      TRUE                                                    ~ NA_character_
    )
  )

# 6. Preserve original sequence_id

db_meta <- rename(db_meta, sequence_id_old = sequence_id)

# 7. Export DB-ready file

fwrite(
  db_meta,
  "/doctorai/chiarba/dataset/stern_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)
