#!/usr/bin/env Rscript

# perez_salvidar_preprocess_for_db.R
# Preprocessing of Perez-Salvidar heavy (IGH) and light (IGK/IGL) chain data for DB.

library(data.table)
library(dplyr)
library(stringr)

# 1. Read IGH AIRR files and tag with SRA_id 

folder_igh  <- "/doctorai/chiarba/add_db/Perez_salvidar/IGH/AIRR"
files_igh   <- list.files(folder_igh, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  base <- tools::file_path_sans_ext(basename(file_path))
  df   <- fread(file_path, colClasses = "character")
  df$file   <- base
  # Example: "result_12345_IGH" -> SRA_id = "12345"
  df$SRA_id <- str_extract(base, "(?<=result_)[0-9]+")
  df
}

db_igh <- rbindlist(lapply(files_igh, read_and_tag), use.names = TRUE, fill = TRUE)

# 2. Read IGK/IGL AIRR files and tag with SRA_id in the same way

folder_igl <- "/doctorai/chiarba/add_db/Perez_salvidar/IGK_IGL/AIRR"
files_igl  <- list.files(folder_igl, pattern = "\\.tsv$", full.names = TRUE)

db_igl <- rbindlist(lapply(files_igl, read_and_tag), use.names = TRUE, fill = TRUE)

# 3. Combine heavy and light chain tables

db <- rbindlist(list(db_igh, db_igl), use.names = TRUE, fill = TRUE)

# 4. Read metadata and prepare minimal set of columns

metadata <- fread("/doctorai/chiarba/add_db/Perez_salvidar/metadata_perez.csv")

metadata <- metadata %>%
  mutate(
    SRA_id = as.character(patient),
    tissue = "blood",
    disease = "RRMS"
  )

meta_subset <- metadata %>%
  select(
    SRA_id,
    tissue,
    disease,
    tretment,      
    BioSample,
    BioProject
  )

# 5. Merge AIRR table with metadata

db$SRA_id <- as.character(db$SRA_id)

db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")

# 6. Remove TCR annotations (both heavy and light, if present)

any(grepl("^TR", db_meta$v_call))
any(grepl("^TR", db_meta$d_call))
any(grepl("^TR", db_meta$j_call))

db_meta <- db_meta %>%
  filter(
    !grepl("^TR", ifelse(is.na(v_call), "", v_call), ignore.case = TRUE),
    !grepl("^TR", ifelse(is.na(d_call), "", d_call), ignore.case = TRUE),
    !grepl("^TR", ifelse(is.na(j_call), "", j_call), ignore.case = TRUE),
    # extra safety: if a 'locus' column exists and starts with TR (TRA/TRB/TRG/TRD), remove that too
    if ("locus" %in% names(.)) !grepl("^TR", ifelse(is.na(locus), "", locus), ignore.case = TRUE) else TRUE
  )

# 7. Create subject_number, study and disease_diagnosis

db_meta <- db_meta %>%
  mutate(
    subject_number = paste0("P_S", SRA_id),
    study          = "Perez_Salvidar"
  ) %>%
  rename(disease_diagnosis = disease)

# 8. Preserve original sequence_id and write DB-ready table

db_meta <- rename(db_meta, sequence_id_old = sequence_id)

fwrite(
  db_meta,
  "/doctorai/chiarba/dataset/PS_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)
