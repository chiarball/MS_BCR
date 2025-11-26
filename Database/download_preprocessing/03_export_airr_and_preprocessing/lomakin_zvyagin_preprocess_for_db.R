#!/usr/bin/env Rscript

# lomakin_zvyagin_preprocess_for_db.R
# The two studies partially overlap (they use some of the same samples and some different ones)
# Preprocessing of Lomakin + Zvyagin AIRR datasets for unified DB.

library(data.table)
library(dplyr)
library(stringr)
library(tools)

## 1. Read AIRR outputs for Zvyagin and Lomakin 

# Zvyagin
path_zvy  <- "/doctorai/chiarba/add_db/Zvyagin/AIRR"
files_zvy <- list.files(path_zvy, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- file_path_sans_ext(basename(file_path))
  df <- fread(file_path, colClasses = "character")  # force all as character
  df$SRA_id <- sra_id
  df
}

zvyagin_sra_db <- rbindlist(lapply(files_zvy, read_and_tag),
                            use.names = TRUE, fill = TRUE)

# Lomakin
path_lom  <- "/doctorai/chiarba/add_db/lomakin/AIRR"
files_lom <- list.files(path_lom, pattern = "\\.tsv$", full.names = TRUE)

lomakin_sra_db <- rbindlist(lapply(files_lom, read_and_tag),
                            use.names = TRUE, fill = TRUE)


## 2. Keep only selected Zvyagin SRA and add metadata 

# These are the Zvyagin SRA_ids you decided to keep (unique or representative)
keep_ids <- c("ERR6501933", "ERR6501936", "ERR6501956", "ERR6501958")

zvy_subset <- zvyagin_sra_db %>%
  filter(SRA_id %in% keep_ids)

# Zvyagin metadata
metadata_zvy <- fread("/doctorai/chiarba/add_db/Zvyagin/metadata_zvy.csv") %>%
  rename(SRA_id = Run)

desired_cols_zvy <- c(
  "SRA_id", "AGE", "Bases", "BioProject", "Broker_name", "cell_type",
  "disease", "disease_staging", "individual", "Organism_part",
  "Sample Name", "Sample_name", "sex", "SRA Study", "disease_duration"
)

metadata_keep_zvy <- metadata_zvy %>%
  filter(SRA_id %in% keep_ids) %>%
  select(any_of(desired_cols_zvy)) %>%
  distinct(SRA_id, .keep_all = TRUE)

zvy_subset_with_meta <- zvy_subset %>%
  left_join(metadata_keep_zvy, by = "SRA_id") %>%
  mutate(study = "zvyagin")


## 3. Add metadata to Lomakin 

metadata_lomakin <- fread("/doctorai/chiarba/add_db/lomakin/lomakin_metadata.csv") %>%
  rename(SRA_id = Run)

desired_cols_loma <- c(
  "SRA_id", "AGE", "cell_type", "common_name", "Consent", "disease",
  "Sample Name", "sex", "SRA Study", "Submitter_Id", "Library Name",
  "disease_duration", "disease_staging", "EDSS", "treatment",
  "BioProject", "BioSample"
)

metadata_lomakin_keep <- metadata_lomakin %>%
  select(any_of(desired_cols_loma)) %>%
  distinct(SRA_id, .keep_all = TRUE)

lomakin_with_meta <- lomakin_sra_db %>%
  left_join(metadata_lomakin_keep, by = "SRA_id")


## 4. Bind Lomakin + Zvyagin in a single table 

lomakin_plus_zvy <- bind_rows(
  zvy_subset_with_meta,
  lomakin_with_meta
) %>%
  mutate(study = dplyr::coalesce(study, "Lomakin"))

db_meta <- lomakin_plus_zvy

# Ensure 'study' has no NA or empty values (safety)
db_meta <- db_meta %>%
  mutate(
    study = if_else(
      is.na(study) | trimws(as.character(study)) == "",
      "Lomakin",
      as.character(study)
    )
  )


## 5. Remove TCR annotations (TR*) 

db_meta <- db_meta %>%
  filter(
    !str_detect(ifelse(is.na(v_call), "", v_call), "^TR"),
    !str_detect(ifelse(is.na(d_call), "", d_call), "^TR"),
    !str_detect(ifelse(is.na(j_call), "", j_call), "^TR")
  )


## 6. Construct subject_number and harmonize clinical fields 

# subject_number rules:
#  - Zvyagin: individual == "MS3" -> "Z3", "MS5" -> "Z5"
#  - Lomakin: from Library Name, "MS<NUMBER>..." -> "Lo<NUMBER>"
db_meta <- db_meta %>%
  mutate(
    .ms_num = str_extract(`Library Name`, "(?<=^MS)\\d+"),
    .lo_from_lib = if_else(!is.na(.ms_num), paste0("Lo", .ms_num), NA_character_),
    subject_number = case_when(
      individual == "MS3" ~ "Z3",
      individual == "MS5" ~ "Z5",
      TRUE                ~ .lo_from_lib
    )
  ) %>%
  select(-.ms_num, -.lo_from_lib)

# treatment: if study is NA/empty (non dovrebbe più succedere), set "no therapy"
db_meta <- db_meta %>%
  mutate(
    treatment = if_else(
      is.na(study) | trimws(as.character(study)) == "",
      "no therapy",
      as.character(treatment)
    )
  )

# disease_diagnosis: all MS
db_meta <- db_meta %>%
  mutate(disease_diagnosis = "MS")

# tissue: rename Organism_part to tissue, then homogenize to "blood"
db_meta <- db_meta %>%
  rename(tissue = Organism_part) %>%
  mutate(tissue = "blood")


## 7. Preserve original sequence_id and write out 

db_meta <- rename(db_meta, sequence_id_old = sequence_id)

fwrite(
  db_meta,
  "/doctorai/chiarba/dataset/lo_zy_for_db.tsv",
  sep = "\t",
  row.names = FALSE
)
