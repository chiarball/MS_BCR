#!/usr/bin/env Rscript

# ramesh_preprocess_for_db.R
# Preprocessing of Ramesh single-cell study for DB.

library(data.table)
library(dplyr)
library(stringr)

# 1. Load AIRR data and metadata 

db <- fread("/doctorai/chiarba/bcr_from_raw/ramesh_final/AIRR_format/ramesh_sra_meta_MSONLY.tsv")

metadata <- fread("/doctorai/chiarba/bcr_from_raw/ramesh/metadata_ramesh.csv")

metadata <- metadata %>%
  rename(SRA_id = Run)

# Select useful columns from metadata
meta_subset <- metadata %>%
  select(
    SRA_id,
    `Library Name`,
    BioProject,
    BioSample
  )

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

# subject_number from subject_id: "Ra" + numeric part
db_meta <- db_meta %>%
  mutate(subject_number = str_c("Ra", str_extract(subject_id, "\\d+")))

db_meta <- db_meta %>%
  mutate(study = "Ramesh")

# (Optional sanity checks)
unique(db_meta$tissue)
unique(db_meta$disease_diagnosis)

# 4. Single-cell cleaning 
# 4.1 Within each SRA_id–cell_id, keep only cells with at most 2 reads

db_meta <- db_meta %>%
  group_by(SRA_id, cell_id) %>%
  mutate(n_per_cell = n()) %>%
  ungroup() %>%
  filter(n_per_cell <= 2) %>%
  select(-n_per_cell)

# 4.2 Within each SRA_id–clone_id, keep only the first occurrence (1 row/clone)

db_meta <- db_meta %>%
  group_by(SRA_id, clone_id) %>%
  mutate(n_cells = n()) %>%
  ungroup() %>%
  distinct(SRA_id, clone_id, .keep_all = TRUE)

# 5. Preserve original sequence_id and write out DB-ready table 

db_meta <- rename(db_meta, sequence_id_old = sequence_id)

fwrite(
  db_meta,
  file = "/doctorai/chiarba/dataset/ramesh_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)
