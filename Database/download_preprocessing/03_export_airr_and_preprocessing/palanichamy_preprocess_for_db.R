#!/usr/bin/env Rscript

# palanichemy_preprocess_for_db.R
# Preprocessing of Palanichemy study for DB.

library(data.table)
library(dplyr)
library(stringr)

# 1. Load AIRR data and metadata

db <- fread("/doctorai/chiarba/bcr_from_raw/pala/fastq_download/AIRR/pala_sra.tsv")

metadata <- fread("/doctorai/chiarba/bcr_from_raw/pala/fastq_download/meta_pala_SRA.csv")

metadata <- metadata %>%
  rename(SRA_id = Run)

# Select useful columns from metadata
meta_subset <- metadata %>%
  select(
    SRA_id,
    sex,
    tissue,
    AGE,
    `Sample Name`,
    isolate,
    BioSample,
    BioProject
  )

# Merge AIRR table with metadata
db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")

# 2. Remove TCR annotations (keep only BCR)

any(grepl("^TR", db$v_call))
any(grepl("^TR", db$d_call))
any(grepl("^TR", db$j_call))

db_meta <- db_meta[
  !(
    grepl("^TR", ifelse(is.na(v_call), "", v_call), ignore.case = TRUE) |
    grepl("^TR", ifelse(is.na(d_call), "", d_call), ignore.case = TRUE) |
    grepl("^TR", ifelse(is.na(j_call), "", j_call), ignore.case = TRUE)
  )
]

# 3. Create subject_number and study columns

# map patient IDs using Sample Name
subject_map <- c(
  "26712" = "Pa1",
  "29612" = "Pa2",
  "14711" = "Pa3",
  "30512" = "Pa4",
  "31012" = "Pa5",
  "34012" = "Pa6",
  "43113" = "Pa7",
  "43213" = "Pa8"
)

db_meta <- db_meta %>%
  mutate(
    subject_number = subject_map[sub("(_.*)$", "", `Sample Name`)]
  )

db_meta <- db_meta %>%
  mutate(study = "Palanichemy")

# 4. Harmonize tissue and disease_diagnosis annotations

# tissue
db_meta <- db_meta %>%
  mutate(tissue = case_when(
    tissue == "cerebrospinal fluid" ~ "CSF",
    tissue == "peripheral blood"    ~ "blood",
    TRUE                            ~ tissue
  ))

# disease diagnosis
diagnosis_map <- c(
  "14711" = "RRMS",
  "26712" = "RRMS",
  "29612" = "RRMS",
  "30512" = "RRMS",
  "31012" = "RRMS",
  "34012" = "PPMS",
  "43113" = "PPMS",
  "43213" = "RRMS"
)

db_meta <- db_meta %>%
  mutate(disease_diagnosis = diagnosis_map[sub("(_.*)$", "", `Sample Name`)])

# 5. Preserve original sequence_id and write out DB-ready table

db_meta <- rename(db_meta, sequence_id_old = sequence_id)

fwrite(
  db_meta,
  file = "/doctorai/chiarba/dataset/pala_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)
