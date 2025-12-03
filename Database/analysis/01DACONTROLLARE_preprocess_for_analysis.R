## Libraries
library(data.table, lib.loc = "/doctorai/chiarba/lib")
library(dplyr,      lib.loc = "/doctorai/chiarba/lib")
library(stringr,    lib.loc = "/doctorai/chiarba/lib")

###############################################################################
## Preprocessing Ryback
###############################################################################

folder_path <- "/doctorai/chiarba/add_db/Ryback/AIRR"
tsv_files   <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df     <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

db <- rbindlist(lapply(tsv_files, read_and_tag), use.names = TRUE, fill = TRUE)

fwrite(db, "/doctorai/chiarba/add_db/Ryback/ryback.tsv",
       sep = "\t", row.names = FALSE)

unique(db$SRA_id)

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
  select(SRA_id, tissue_type, `Sample Name`, `Library Name`, sex, agegroup)

db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")

fwrite(db_meta, "/doctorai/chiarba/add_db/Ryback/ryback_metadata.tsv",
       sep = "\t", row.names = FALSE)

db_meta_sel <- db_meta %>%
  select(
    SRA_id,
    sequence,
    junction,
    sequence_alignment,
    germline_alignment,
    cdr3,
    cdr3_aa,
    v_call,
    d_call,
    j_call,
    junction_aa,
    junction_length,
    `Sample Name`,
    `Library Name`
  ) %>%
  mutate(
    v_gene = sub("\\*.*", "", v_call),
    d_gene = sub("\\*.*", "", d_call),
    j_gene = sub("\\*.*", "", j_call)
  )

db_meta_sel <- subset(db_meta_sel, grepl("^IGH", v_gene))

db_meta_sel <- db_meta_sel %>%
  filter(!is.na(junction_aa) & junction_aa != "") %>%
  filter(!grepl("\\*", junction_aa)) %>%
  filter(!grepl("_",  junction_aa))

db_meta_sel <- db_meta_sel[grepl("^C.*W$", db_meta_sel$junction_aa), ]

any(grepl("^TR", db_meta_sel$v_call))
any(grepl("^TR", db_meta_sel$d_gene))
any(grepl("^TR", db_meta_sel$j_call))

db_meta_sel <- db_meta_sel %>%
  mutate(
    v_cdr3_j = paste(v_gene, junction_aa, j_gene, sep = "_"),
    subject_number = str_extract(`Library Name`, "(?<=SBL0)\\d+") %>%
      paste0("Ry", .),
    study = "Ryback",
    tissue = "blood",
    tissue_type = case_when(
      tissue == "blood" ~ "extra",
      tissue == "CSF"   ~ "intra"
    ),
    junction_aa_length = as.numeric(junction_length) / 3,
    disease_diagnosis  = "MS"
  )

fwrite(db_meta_sel,
       "/doctorai/chiarba/add_db/Ryback/ryback_for_my_analysis.tsv",
       sep = "\t", row.names = FALSE)

###############################################################################
## Preprocessing Perez_Salvidar (IGH)
###############################################################################

folder_path <- "/doctorai/chiarba/add_db/Perez_salvidar/IGH/AIRR"
tsv_files   <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df     <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

db <- rbindlist(lapply(tsv_files, read_and_tag), use.names = TRUE, fill = TRUE)

fwrite(db,
       "/doctorai/chiarba/add_db/Perez_salvidar/IGH/Perez_salvidar_IGH.tsv",
       sep = "\t", row.names = FALSE)

db <- db %>%
  rename(file = SRA_id) %>%
  mutate(SRA_id = str_extract(file, "(?<=result_)[0-9]+"))

unique(db$SRA_id)

metadata <- fread("/doctorai/chiarba/add_db/Perez_salvidar/metadata_perez.csv") %>%
  mutate(
    SRA_id = patient,
    tissue = "blood",
    disease = "RRMS"
  )

meta_subset <- metadata %>%
  select(SRA_id, tissue, disease, tretment)

db$SRA_id         <- as.character(db$SRA_id)
meta_subset$SRA_id <- as.character(meta_subset$SRA_id)

db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")

fwrite(db_meta,
       "/doctorai/chiarba/add_db/Perez_salvidar/IGH/perez_sal_IGH_metadata.tsv",
       sep = "\t", row.names = FALSE)

db_meta_sel <- db_meta %>%
  select(
    SRA_id,
    sequence,
    junction,
    sequence_alignment,
    germline_alignment,
    cdr3,
    cdr3_aa,
    v_call,
    d_call,
    j_call,
    junction_aa,
    junction_length,
    file,
    tissue,
    disease,
    tretment
  ) %>%
  mutate(
    v_gene = sub("\\*.*", "", v_call),
    d_gene = sub("\\*.*", "", d_call),
    j_gene = sub("\\*.*", "", j_call)
  )

db_meta_sel <- subset(db_meta_sel, grepl("^IGH", v_gene))

db_meta_sel <- db_meta_sel %>%
  filter(!is.na(junction_aa) & junction_aa != "") %>%
  filter(!grepl("\\*", junction_aa)) %>%
  filter(!grepl("_",  junction_aa))

db_meta_sel <- db_meta_sel[grepl("^C.*W$", db_meta_sel$junction_aa), ]

any(grepl("^TR", db_meta_sel$v_call))
any(grepl("^TR", db_meta_sel$d_gene))
any(grepl("^TR", db_meta_sel$j_call))

db_meta_sel <- db_meta_sel %>%
  mutate(
    v_cdr3_j = paste(v_gene, junction_aa, j_gene, sep = "_"),
    subject_number = paste0("P_S", SRA_id),
    study = "Perez_Salvidar",
    tissue_type = case_when(
      tissue == "blood" ~ "extra",
      tissue == "CSF"   ~ "intra"
    ),
    junction_aa_length = as.numeric(junction_length) / 3
  ) %>%
  rename(disease_diagnosis = disease)

unique(db_meta_sel$tretment)

db_meta_sel_naive <- db_meta_sel %>%
  filter(tretment == "naive")

unique(db_meta_sel_naive$subject_number)

fwrite(db_meta_sel_naive,
       "/doctorai/chiarba/add_db/Perez_salvidar/IGH/perez_naive_for_my_analysis.tsv",
       sep = "\t", row.names = FALSE)

###############################################################################
## Preprocessing Stern
###############################################################################

folder_path <- "/doctorai/chiarba/bcr_from_raw/stern/fastq_download/AIRR"
tsv_files   <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df     <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

stern_sra_db <- rbindlist(lapply(tsv_files, read_and_tag),
                          use.names = TRUE, fill = TRUE)

fwrite(stern_sra_db,
       "/doctorai/chiarba/bcr_from_raw/stern/fastq_download/AIRR/stern_sra_db.tsv",
       sep = "\t", row.names = FALSE)

metadata_stern <- fread("/doctorai/chiarba/bcr_from_raw/stern/meta_stern.csv") %>%
  filter(isolate != "M1") %>%
  rename(SRA_id = Run)

meta_subset <- metadata_stern %>%
  select(SRA_id, sex, tissue, AGE, `Sample Name`, disease_stage)

stern_sra_db_enriched <- stern_sra_db %>%
  left_join(meta_subset, by = "SRA_id")

fwrite(stern_sra_db_enriched,
       "/doctorai/chiarba/bcr_from_raw/stern/fastq_download/AIRR/stern_sra_db.tsv",
       sep = "\t", row.names = FALSE)

stern_sra_selected <- stern_sra_db_enriched %>%
  select(
    SRA_id,
    sequence,
    junction,
    sequence_alignment,
    germline_alignment,
    cdr3,
    cdr3_aa,
    v_call,
    d_call,
    j_call,
    junction_aa,
    tissue,
    junction_length,
    `Sample Name`,
    disease_stage
  ) %>%
  mutate(
    v_gene = sub("\\*.*", "", v_call),
    d_gene = sub("\\*.*", "", d_call),
    j_gene = sub("\\*.*", "", j_call)
  )

stern_sra_selected <- subset(stern_sra_selected, grepl("^IGH", v_gene))

stern_sra_selected <- stern_sra_selected %>%
  filter(!is.na(junction_aa) & junction_aa != "") %>%
  filter(!grepl("\\*", junction_aa)) %>%
  filter(!grepl("_",  junction_aa))

stern_sra_selected <- stern_sra_selected[grepl("^C.*W$", stern_sra_selected$junction_aa), ]

any(grepl("^TR", stern_sra_selected$d_gene))
any(grepl("^TR", stern_sra_selected$v_call))
any(grepl("^TR", stern_sra_selected$j_call))

stern_sra_selected <- stern_sra_selected %>%
  mutate(
    v_cdr3_j = paste(v_gene, junction_aa, j_gene, sep = "_"),
    subject_number = str_extract(`Sample Name`, "^M\\d+"),
    study   = "Stern",
    tissue  = case_when(
      tissue == "Cervical lymph node" ~ "CLN",
      TRUE                            ~ tissue
    ),
    tissue_type = case_when(
      str_detect(`Sample Name`, "CLN") ~ "extra",
      TRUE                             ~ "intra"
    ),
    junction_aa_length = as.numeric(junction_length) / 3,
    disease_diagnosis = case_when(
      disease_stage == "Primary progressive"                 ~ "PP",
      disease_stage %in% c("Chronic progressive", "Chronic") ~ "Chronic",
      disease_stage == "Secondary progressive"               ~ "SP",
      TRUE                                                   ~ NA_character_
    )
  )

any(is.na(stern_sra_selected$disease_diagnosis) |
      stern_sra_selected$disease_diagnosis == "")

stern_sra_selected$subject_number <- sub("^M", "S", stern_sra_selected$subject_number)

fwrite(stern_sra_selected,
       "/doctorai/chiarba/bcr_from_raw/stern/fastq_download/AIRR/stern_cleaned_selected_igh.tsv",
       sep = "\t", row.names = FALSE)

###############################################################################
## Preprocessing Agrafiotis
###############################################################################

path  <- "/doctorai/chiarba/add_db/agrafiotis/AIRR_format"
files <- list.files(path, pattern = "\\.tsv$", full.names = TRUE)

agrafiotis <- rbindlist(
  lapply(files, function(f) {
    dt <- fread(f, sep = "\t", na.strings = c("", "NA"))
    dt[, `__source_file` := basename(f)]
    dt
  }),
  fill = TRUE
)

fwrite(agrafiotis,
       "/doctorai/chiarba/add_db/agrafiotis/agrafiotis_all.tsv",
       sep = "\t", row.names = FALSE)

metadata_arg <- fread("/doctorai/chiarba/add_db/agrafiotis/metadataagrafiotis.csv") %>%
  rename(
    SRA_id   = Run,
    sample_id = `Sample Name`
  )

agrafiotis$SRR <- sub(".*_(SRR[0-9]+)\\.contigs\\.tsv$", "\\1",
                      agrafiotis$`__source_file`)

agrafiotis <- agrafiotis %>%
  mutate(SRA_id = str_extract(`__source_file`, "SRR\\d+"))

meta_subset <- metadata_arg %>%
  select(SRA_id, sex, tissue, AGE, disease_stage, sample_id)

agr_sra_meta <- agrafiotis %>%
  left_join(meta_subset, by = "SRA_id")

fwrite(agr_sra_meta,
       "/doctorai/chiarba/add_db/agrafiotis/agrafiotis_all_meta.tsv",
       sep = "\t", row.names = FALSE)

agr_sra_meta <- fread("/doctorai/chiarba/add_db/agrafiotis/agrafiotis_all_meta.tsv")

agr_sra_meta <- agr_sra_meta %>%
  group_by(SRA_id, cell_id) %>%
  mutate(n_per_cell = n()) %>%
  ungroup() %>%
  filter(n_per_cell <= 2) %>%
  select(-n_per_cell)

fwrite(agr_sra_meta,
       "/doctorai/chiarba/add_db/agrafiotis/agr_cellid.tsv",
       sep = "\t", row.names = FALSE)

clone_counts <- agr_sra_meta %>%
  group_by(SRA_id, clone_id) %>%
  summarise(
    n_occurrences = n(),
    .groups = "drop"
  )

## View(clone_counts)  # optional interactive check

prova <- agr_sra_meta %>%
  filter(SRA_id == "SRR23130223", clone_id == 0) %>%
  summarise(n_distinct_cdr3 = n_distinct(sequence))

agr_collapsed <- agr_sra_meta %>%
  group_by(SRA_id, clone_id) %>%
  mutate(n_cells = n()) %>%
  ungroup() %>%
  distinct(SRA_id, clone_id, .keep_all = TRUE)

fwrite(agr_collapsed,
       "/doctorai/chiarba/add_db/agrafiotis/agr_cellid_cloneid.tsv",
       sep = "\t", row.names = FALSE)

agr_sra_meta <- fread("/doctorai/chiarba/add_db/agrafiotis/agr_cellid_cloneid.tsv")

agr_sra_meta <- agr_sra_meta %>%
  select(
    SRA_id,
    sequence,
    junction,
    sequence_alignment,
    germline_alignment,
    cdr3,
    cdr3_aa,
    v_call,
    d_call,
    j_call,
    junction_aa,
    tissue,
    junction_length,
    sample_id,
    disease_stage,
    clone_id,
    sequence_id,
    cell_ids,
    cell_id
  ) %>%
  mutate(
    v_gene = sub("\\*.*", "", v_call),
    d_gene = sub("\\*.*", "", d_call),
    j_gene = sub("\\*.*", "", j_call)
  ) %>%
  filter(!is.na(junction_aa) & junction_aa != "") %>%
  filter(!grepl("\\*", junction_aa)) %>%
  filter(!grepl("_",  junction_aa))

any(grepl("^TR", agr_sra_meta$d_gene))
any(grepl("^TR", agr_sra_meta$v_call))
any(grepl("^TR", agr_sra_meta$j_call))

agr_sra_meta <- agr_sra_meta %>%
  filter(
    !str_detect(d_gene, "^TR"),
    !str_detect(v_call, "^TR"),
    !str_detect(j_call, "^TR")
  )

agr_sra_meta <- agr_sra_meta %>%
  mutate(
    subject_number = paste0("A", substr(sample_id, 1, 1),
                            sub(".*-(\\d+)$", "\\1", sample_id)),
    study = "Agrafiotis",
    tissue_type = case_when(
      str_detect(tissue, regex("CSF", ignore_case = TRUE))   ~ "intra",
      str_detect(tissue, regex("blood", ignore_case = TRUE)) ~ "extra",
      TRUE                                                   ~ NA_character_
    ),
    junction_aa_length = as.numeric(junction_length) / 3
  )

agr_sra_meta <- subset(agr_sra_meta, grepl("^IGH", v_gene))
agr_sra_meta <- agr_sra_meta[grepl("^C.*W$", agr_sra_meta$junction_aa), ]

agr_sra_meta <- agr_sra_meta %>%
  mutate(
    v_cdr3_j      = paste(v_gene, junction_aa, j_gene, sep = "_")
  ) %>%
  rename(
    mixcr_clone_id     = clone_id,
    disease_diagnosis  = disease_stage
  )

fwrite(agr_sra_meta,
       "/doctorai/chiarba/add_db/agrafiotis/agr_colse_for_my_analysis.tsv",
       sep = "\t", row.names = FALSE)

###############################################################################
## Preprocessing Greenfield
###############################################################################

folder_path <- "/doctorai/chiarba/bcr_from_raw/green/AIRR"
tsv_files   <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df     <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

Green_db <- rbindlist(lapply(tsv_files, read_and_tag),
                      use.names = TRUE, fill = TRUE)

fwrite(Green_db,
       "/doctorai/chiarba/bcr_from_raw/green/AIRR/Green_all.tsv",
       sep = "\t", row.names = FALSE)

metadata_green <- fread("/doctorai/chiarba/bcr_from_raw/green/AIRR/Metadata_green_SraRunTable.csv") %>%
  rename(SRA_id = Run)

meta_subset <- metadata_green %>%
  select(SRA_id, sex, tissue, age, disease, `Library Name`)

Green_db_enriched2 <- Green_db %>%
  left_join(meta_subset, by = "SRA_id")

fwrite(Green_db_enriched2,
       "/doctorai/chiarba/bcr_from_raw/green/AIRR/Green_all_meta.tsv",
       sep = "\t", row.names = FALSE)

Green_db_selected <- Green_db_enriched2 %>%
  select(
    SRA_id,
    sequence,
    junction,
    sequence_alignment,
    germline_alignment,
    cdr3,
    cdr3_aa,
    v_call,
    d_call,
    j_call,
    junction_aa,
    disease,
    tissue,
    junction_length,
    `Library Name`
  ) %>%
  mutate(
    v_gene = sub("\\*.*", "", v_call),
    d_gene = sub("\\*.*", "", d_call),
    j_gene = sub("\\*.*", "", j_call)
  )

Green_db_selected <- subset(Green_db_selected, grepl("^IGH", v_gene))

Green_db_selected <- Green_db_selected %>%
  filter(!is.na(junction_aa) & junction_aa != "") %>%
  filter(!grepl("\\*", junction_aa)) %>%
  filter(!grepl("_",  junction_aa))

Green_db_selected <- Green_db_selected[grepl("^C.*W$", Green_db_selected$junction_aa), ]

Green_db_selected <- Green_db_selected %>%
  mutate(
    v_cdr3_j = paste(v_gene, junction_aa, j_gene, sep = "_")
  )

any(grepl("^TR", Green_db_selected$v_gene))
any(grepl("^TR", Green_db_selected$d_gene))
any(grepl("^TR", Green_db_selected$j_gene))

Green_db_selected <- Green_db_selected %>%
  rename(subject_number = `Library Name`) %>%
  mutate(
    subject_number = str_extract(subject_number, "^Pt\\d+"),
    subject_number = str_replace(subject_number, "^Pt", "G"),
    study          = "Greenfield",
    tissue = case_when(
      tissue == "cerebrospinal fluid" ~ "CSF",
      tissue == "peripheral blood"    ~ "blood",
      TRUE                            ~ tissue
    ),
    disease_diagnosis = case_when(
      disease == "multiple sclerosis" ~ "MS",
      TRUE                            ~ disease
    ),
    tissue_type = case_when(
      tissue == "blood" ~ "extra",
      tissue == "CSF"   ~ "intra"
    ),
    junction_aa_length = as.numeric(junction_length) / 3
  ) %>%
  select(-junction_length)

fwrite(Green_db_selected,
       "/doctorai/chiarba/bcr_from_raw/green/AIRR/Green_cleaned_selectedcol_igh.tsv",
       sep = "\t", row.names = FALSE)

###############################################################################
## Preprocessing Palanichemy
###############################################################################

folder_path <- "/doctorai/chiarba/bcr_from_raw/pala/fastq_download/AIRR"
tsv_files   <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df     <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

pala_sra_db <- rbindlist(lapply(tsv_files, read_and_tag),
                         use.names = TRUE, fill = TRUE)

fwrite(pala_sra_db,
       "/doctorai/chiarba/bcr_from_raw/pala/fastq_download/AIRR/pala_sra.tsv",
       sep = "\t", row.names = FALSE)

metadata_pala <- fread("/doctorai/chiarba/bcr_from_raw/pala/fastq_download/meta_pala_SRA.csv") %>%
  rename(SRA_id = Run)

meta_subset <- metadata_pala %>%
  select(SRA_id, sex, tissue, AGE, `Sample Name`, isolate, BioSample, BioProject)

pala_sra_db_enriched <- pala_sra_db %>%
  left_join(meta_subset, by = "SRA_id")

fwrite(pala_sra_db_enriched,
       "/doctorai/chiarba/bcr_from_raw/pala/fastq_download/pala_sra_meta.tsv",
       sep = "\t", row.names = FALSE)

pala_sra_selected <- pala_sra_db_enriched %>%
  select(
    SRA_id,
    sequence,
    junction,
    sequence_alignment,
    germline_alignment,
    cdr3,
    cdr3_aa,
    v_call,
    d_call,
    j_call,
    junction_aa,
    tissue,
    junction_length,
    `Sample Name`
  ) %>%
  mutate(
    v_gene = sub("\\*.*", "", v_call),
    d_gene = sub("\\*.*", "", d_call),
    j_gene = sub("\\*.*", "", j_call)
  )

pala_sra_selected <- subset(pala_sra_selected, grepl("^IGH", v_gene))

pala_sra_selected <- pala_sra_selected %>%
  filter(!is.na(junction_aa) & junction_aa != "") %>%
  filter(!grepl("\\*", junction_aa)) %>%
  filter(!grepl("_",  junction_aa))

pala_sra_selected <- pala_sra_selected[grepl("^C.*W$", pala_sra_selected$junction_aa), ]

any(grepl("^TR", pala_sra_selected$v_call))
any(grepl("^TR", pala_sra_selected$d_gene))
any(grepl("^TR", pala_sra_selected$j_call))

pala_sra_selected <- pala_sra_selected %>%
  mutate(
    v_cdr3_j = paste(v_gene, junction_aa, j_gene, sep = "_")
  )

subject_map <- c(
  "26712" = "P1",
  "29612" = "P2",
  "14711" = "P3",
  "30512" = "P4",
  "31012" = "P5",
  "34012" = "P6",
  "43113" = "P7",
  "43213" = "P8"
)

pala_sra_selected <- pala_sra_selected %>%
  mutate(
    subject_number = subject_map[sub("(_.*)$", "", `Sample Name`)],
    study  = "Palanichemy",
    tissue = case_when(
      tissue == "cerebrospinal fluid" ~ "CSF",
      tissue == "peripheral blood"    ~ "blood",
      TRUE                            ~ tissue
    ),
    tissue_type = case_when(
      tissue == "blood" ~ "extra",
      tissue == "CSF"   ~ "intra"
    ),
    junction_aa_length = as.numeric(junction_length) / 3
  )

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

pala_sra_selected <- pala_sra_selected %>%
  mutate(
    disease_diagnosis = diagnosis_map[sub("(_.*)$", "", `Sample Name`)]
  )

fwrite(pala_sra_selected,
       "/doctorai/chiarba/bcr_from_raw/pala/pala_cleaned_selectedcol_igh.tsv",
       sep = "\t", row.names = FALSE)

###############################################################################
## Preprocessing Laurent
###############################################################################

folder_path <- "/doctorai/chiarba/add_db/Laurent/AIRR"
tsv_files   <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df     <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

db <- rbindlist(lapply(tsv_files, read_and_tag),
                use.names = TRUE, fill = TRUE)

fwrite(db,
       "/doctorai/chiarba/add_db/Laurent/laurent.tsv",
       sep = "\t", row.names = FALSE)

unique(db$SRA_id)

metadata <- fread("/doctorai/chiarba/add_db/Laurent/metadata_lau.csv") %>%
  rename(SRA_id = Run)

meta_subset <- metadata %>%
  select(
    SRA_id, sex, tissue, AGE, `Sample Name`, Patient_ID, receptor_type,
    time_point, treatment, weeks_since_baseline
  )

db_meta <- db %>%
  left_join(meta_subset, by = "SRA_id")

fwrite(db_meta,
       "/doctorai/chiarba/add_db/Laurent/laurent_metadata.tsv",
       sep = "\t", row.names = FALSE)

db_meta_sel <- db_meta %>%
  select(
    SRA_id,
    sequence,
    junction,
    sequence_alignment,
    germline_alignment,
    cdr3,
    cdr3_aa,
    v_call,
    d_call,
    j_call,
    junction_aa,
    tissue,
    junction_length,
    `Sample Name`,
    Patient_ID,
    receptor_type,
    time_point,
    treatment,
    weeks_since_baseline
  ) %>%
  mutate(
    v_gene = sub("\\*.*", "", v_call),
    d_gene = sub("\\*.*", "", d_call),
    j_gene = sub("\\*.*", "", j_call)
  )

db_meta_sel <- subset(db_meta_sel, grepl("^IGH", v_gene))

db_meta_sel <- db_meta_sel %>%
  filter(!is.na(junction_aa) & junction_aa != "") %>%
  filter(!grepl("\\*", junction_aa)) %>%
  filter(!grepl("_",  junction_aa))

db_meta_sel <- db_meta_sel[grepl("^C.*W$", db_meta_sel$junction_aa), ]

any(grepl("^TR", db_meta_sel$v_gene))
any(grepl("^TR", db_meta_sel$d_gene))
any(grepl("^TR", db_meta_sel$j_gene))

db_meta_sel <- db_meta_sel %>%
  mutate(
    v_cdr3_j = paste(v_gene, junction_aa, j_gene, sep = "_"),
    subject_number = paste0("Lau", sub("MS-(\\d+)", "\\1", Patient_ID)),
    study  = "Laurent",
    tissue = "blood",
    tissue_type = case_when(
      tissue == "blood" ~ "extra",
      tissue == "CSF"   ~ "intra"
    ),
    junction_aa_length = as.numeric(junction_length) / 3,
    disease_diagnosis  = "MS"
  ) %>%
  filter(weeks_since_baseline == "0")

unique(db_meta_sel$receptor_type)

fwrite(db_meta_sel,
       "/doctorai/chiarba/add_db/Laurent/laurent_naive_IgG_IgM_for_my_analysis.tsv",
       sep = "\t", row.names = FALSE)

###############################################################################
## Preprocessing Ramesh
###############################################################################

folder_path <- "/doctorai/chiarba/bcr_from_raw/ramesh_final/AIRR_format"
tsv_files   <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df     <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

ramesh_raw <- rbindlist(lapply(tsv_files, read_and_tag),
                        use.names = TRUE, fill = TRUE)

unique(ramesh_raw$SRA_id)

ramesh_raw$SRA_id <- str_extract(ramesh_raw$SRA_id, "SRR\\d+")

unique(ramesh_raw$SRA_id)

fwrite(ramesh_raw,
       "/doctorai/chiarba/bcr_from_raw/ramesh_final/AIRR_format/ramesh_sra_db.tsv",
       sep = "\t", row.names = FALSE)

metadata_ramesh <- fread("/doctorai/chiarba/bcr_from_raw/ramesh/metadata_ramesh.csv") %>%
  rename(
    SRA_id   = Run,
    sample_id = `Sample Name`
  ) %>%
  mutate(subject_id = str_extract(sample_id, "^[0-9]+")) %>%
  mutate(
    disease_diagnosis = case_when(
      subject_id %in% c("2", "3", "4", "5", "21", "32") ~ "other",
      subject_id %in% c("31", "22")                     ~ "CIS",
      TRUE                                              ~ "RRMS"
    )
  )

fwrite(metadata_ramesh,
       "/doctorai/chiarba/bcr_from_raw/ramesh/meta_ramesh_DD.tsv",
       sep = "\t", row.names = FALSE)

meta_subset <- metadata_ramesh %>%
  select(SRA_id, sex, tissue, AGE, disease_diagnosis, subject_id)

ramesh_sra_meta <- ramesh_raw %>%
  left_join(meta_subset, by = "SRA_id")

fwrite(ramesh_sra_meta,
       "/doctorai/chiarba/bcr_from_raw/ramesh_final/AIRR_format/ramesh_sra_meta.tsv",
       sep = "\t", row.names = FALSE)

ramesh_sra_meta <- ramesh_sra_meta %>%
  filter(disease_diagnosis %in% c("RRMS", "CIS"))

unique(ramesh_sra_meta$disease_diagnosis)

fwrite(ramesh_sra_meta,
       "/doctorai/chiarba/bcr_from_raw/ramesh_final/AIRR_format/ramesh_sra_meta_MSONLY.tsv",
       sep = "\t", row.names = FALSE)

ramesh_sra_meta <- ramesh_sra_meta %>%
  select(
    SRA_id,
    sequence,
    junction,
    sequence_alignment,
    germline_alignment,
    cdr3,
    cdr3_aa,
    v_call,
    d_call,
    j_call,
    junction_aa,
    tissue,
    junction_length,
    subject_id,
    disease_diagnosis,
    clone_id,
    sequence_id,
    cell_ids,
    cell_id
  ) %>%
  mutate(
    v_gene = sub("\\*.*", "", v_call),
    d_gene = sub("\\*.*", "", d_call),
    j_gene = sub("\\*.*", "", j_call)
  ) %>%
  filter(!is.na(junction_aa) & junction_aa != "") %>%
  filter(!grepl("\\*", junction_aa)) %>%
  filter(!grepl("_",  junction_aa))

any(grepl("^TR", ramesh_sra_meta$d_gene))
any(grepl("^TR", ramesh_sra_meta$v_call))
any(grepl("^TR", ramesh_sra_meta$j_call))

ramesh_sra_meta <- ramesh_sra_meta %>%
  filter(
    !str_detect(d_gene, "^TR"),
    !str_detect(v_call, "^TR"),
    !str_detect(j_call, "^TR")
  ) %>%
  mutate(
    subject_number = str_c("R", str_extract(subject_id, "\\d+")),
    study = "Ramesh",
    tissue_type = case_when(
      str_detect(tissue, regex("CSF", ignore_case = TRUE))   ~ "intra",
      str_detect(tissue, regex("blood", ignore_case = TRUE)) ~ "extra",
      TRUE                                                   ~ NA_character_
    ),
    junction_aa_length = as.numeric(junction_length) / 3
  )

fwrite(ramesh_sra_meta,
       "/doctorai/chiarba/bcr_from_raw/ramesh_final/AIRR_format/ramesh_selected.tsv",
       sep = "\t", row.names = FALSE)

ramesh_sra_meta <- ramesh_sra_meta %>%
  group_by(SRA_id, cell_id) %>%
  mutate(n_per_cell = n()) %>%
  ungroup() %>%
  filter(n_per_cell <= 2) %>%
  select(-n_per_cell)

fwrite(ramesh_sra_meta,
       "/doctorai/chiarba/bcr_from_raw/ramesh_final/AIRR_format/ramesh_selected_cleaned_cellid.tsv",
       sep = "\t", row.names = FALSE)

clone_counts <- ramesh_sra_meta %>%
  group_by(SRA_id, clone_id) %>%
  summarise(
    n_occurrences = n(),
    .groups = "drop"
  )

## View(clone_counts)  # optional interactive check

prova <- ramesh_sra_meta %>%
  filter(SRA_id == "SRR12483435", clone_id == 100) %>%
  summarise(n_distinct_cdr3 = n_distinct(sequence))

ramesh_collapsed <- ramesh_sra_meta %>%
  group_by(SRA_id, clone_id) %>%
  mutate(n_cells = n()) %>%
  ungroup() %>%
  distinct(SRA_id, clone_id, .keep_all = TRUE)

fwrite(ramesh_collapsed,
       "/doctorai/chiarba/bcr_from_raw/ramesh_final/AIRR_format/ramesh_selected_cleaned_cellid_cloneid.tsv",
       sep = "\t", row.names = FALSE)

ramesh_collapsed <- subset(ramesh_collapsed, grepl("^IGH", v_gene))
ramesh_collapsed <- ramesh_collapsed[grepl("^C.*W$", ramesh_collapsed$junction_aa), ]

ramesh_collapsed <- ramesh_collapsed %>%
  mutate(
    v_cdr3_j = paste(v_gene, junction_aa, j_gene, sep = "_")
  ) %>%
  rename(mixcr_clone_id = clone_id)

fwrite(ramesh_collapsed,
       "/doctorai/chiarba/bcr_from_raw/ramesh_final/AIRR_format/ramesh_final_preprocessing.tsv",
       sep = "\t", row.names = FALSE)

###############################################################################
## Preprocessing Zvyagin and Lomakin (merge + selection for analysis)
###############################################################################

folder_path <- "/doctorai/chiarba/add_db/Zvyagin/AIRR"
tsv_files   <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df     <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

zvyagin_sra_db <- rbindlist(lapply(tsv_files, read_and_tag),
                            use.names = TRUE, fill = TRUE)

unique(zvyagin_sra_db$SRA_id)

fwrite(zvyagin_sra_db,
       "/doctorai/chiarba/add_db/Zvyagin/zvyagin_airr_output_only_merged.tsv",
       sep = "\t", row.names = FALSE)

folder_path <- "/doctorai/chiarba/add_db/lomakin/AIRR"
tsv_files   <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df     <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

lomakin_sra_db <- rbindlist(lapply(tsv_files, read_and_tag),
                            use.names = TRUE, fill = TRUE)

fwrite(lomakin_sra_db,
       "/doctorai/chiarba/add_db/lomakin/lomakin_airr_output_only_merged.tsv",
       sep = "\t", row.names = FALSE)

subdbLomakin <- subset(lomakin_sra_db, SRA_id == "ERR4588280")
subdbzvy     <- subset(zvyagin_sra_db, SRA_id == "ERR6501961")

seq_zvy     <- subdbzvy$sequence
seq_lomakin <- subdbLomakin$sequence

n_common    <- length(intersect(seq_zvy, seq_lomakin))
n_only_zvy  <- length(setdiff(seq_zvy, seq_lomakin))
n_only_loma <- length(setdiff(seq_lomakin, seq_zvy))

cat("Common sequences:            ", n_common,    "\n")
cat("Sequences only in subdbzvy:  ", n_only_zvy,  "\n")
cat("Sequences only in subdbLoma: ", n_only_loma, "\n")

lomakin_filtered <- lomakin_sra_db %>%
  filter(SRA_id != "ERR4588287")

zvyagin_filtered <- zvyagin_sra_db %>%
  filter(!SRA_id %in% c("ERR6501933", "ERR6501936", "ERR6501956", "ERR6501958"))

seq_loma <- lomakin_filtered$sequence
seq_zvy  <- zvyagin_filtered$sequence

n_common    <- length(intersect(seq_loma, seq_zvy))
n_only_loma <- length(setdiff(seq_loma, seq_zvy))
n_only_zvy  <- length(setdiff(seq_zvy,  seq_loma))

cat("Number of common sequences:          ", n_common,    "\n")
cat("Number of sequences only in lomakin: ", n_only_loma, "\n")
cat("Number of sequences only in zvyagin: ", n_only_zvy,  "\n")

lomakin_filtered <- subset(lomakin_filtered, !is.na(sequence))
zvyagin_filtered <- subset(zvyagin_filtered, !is.na(sequence))

seq_only_loma <- setdiff(lomakin_filtered$sequence, zvyagin_filtered$sequence)
seq_only_zvy  <- setdiff(zvyagin_filtered$sequence,  lomakin_filtered$sequence)

rows_loma_only <- unique(lomakin_filtered[lomakin_filtered$sequence %in% seq_only_loma,
                                          c("SRA_id", "sequence")])
rows_zvy_only  <- unique(zvyagin_filtered[zvyagin_filtered$sequence %in% seq_only_zvy,
                                          c("SRA_id", "sequence")])

cat("Lomakin-only sequences (first 20 rows):\n")
print(head(rows_loma_only, 20))

cat("\nZvyagin-only sequences (first 20 rows):\n")
print(head(rows_zvy_only, 20))

SRA_loma_only <- sort(unique(rows_loma_only$SRA_id))
SRA_zvy_only  <- sort(unique(rows_zvy_only$SRA_id))

cat("\nDistinct SRA_id with Lomakin-only sequences (n=", length(SRA_loma_only), "):\n", sep = "")
print(SRA_loma_only)

cat("\nDistinct SRA_id with Zvyagin-only sequences (n=", length(SRA_zvy_only), "):\n", sep = "")
print(SRA_zvy_only)

shared_seq  <- intersect(lomakin_filtered$sequence, zvyagin_filtered$sequence)
unique_loma <- setdiff(lomakin_filtered$sequence, zvyagin_filtered$sequence)
unique_zvy  <- setdiff(zvyagin_filtered$sequence,  lomakin_filtered$sequence)

loma_pairs <- lomakin_filtered %>%
  distinct(SRA_id, sequence)
zvy_pairs  <- zvyagin_filtered %>%
  distinct(SRA_id, sequence)

loma_stats <- loma_pairs %>%
  mutate(
    is_unique = sequence %in% unique_loma,
    is_shared = sequence %in% shared_seq
  ) %>%
  group_by(SRA_id) %>%
  summarise(
    n_seq     = n(),
    n_unique  = sum(is_unique),
    n_shared  = sum(is_shared),
    all_unique = n_unique == n_seq,
    .groups   = "drop"
  ) %>%
  arrange(desc(n_unique))

zvy_stats <- zvy_pairs %>%
  mutate(
    is_unique = sequence %in% unique_zvy,
    is_shared = sequence %in% shared_seq
  ) %>%
  group_by(SRA_id) %>%
  summarise(
    n_seq     = n(),
    n_unique  = sum(is_unique),
    n_shared  = sum(is_shared),
    all_unique = n_unique == n_seq,
    .groups   = "drop"
  ) %>%
  arrange(desc(n_unique))

SRA_loma_any_unique <- loma_stats %>%
  filter(n_unique > 0) %>%
  pull(SRA_id)
SRA_zvy_any_unique  <- zvy_stats %>%
  filter(n_unique > 0) %>%
  pull(SRA_id)

SRA_loma_all_unique <- loma_stats %>%
  filter(all_unique) %>%
  pull(SRA_id)
SRA_zvy_all_unique <- zvy_stats %>%
  filter(all_unique) %>%
  pull(SRA_id)

cat("Lomakin: SRA_id with any unique seq:\n"); print(SRA_loma_any_unique)
cat("\nLomakin: SRA_id with ALL sequences unique:\n"); print(SRA_loma_all_unique)

cat("\nZvyagin: SRA_id with any unique seq:\n"); print(SRA_zvy_any_unique)
cat("\nZvyagin: SRA_id with ALL sequences unique:\n"); print(SRA_zvy_all_unique)

lomakin_pairs <- lomakin_filtered %>%
  filter(!is.na(sequence)) %>%
  distinct(SRA_id, sequence)

zvy_pairs <- zvyagin_filtered %>%
  filter(!is.na(sequence)) %>%
  distinct(SRA_id, sequence)

shared_seq <- intersect(lomakin_pairs$sequence, zvy_pairs$sequence)

summarise_per_sra <- function(pairs_df) {
  pairs_df %>%
    mutate(
      is_shared = sequence %in% shared_seq,
      is_unique = !is_shared
    ) %>%
    group_by(SRA_id) %>%
    summarise(
      n_total  = n(),
      n_shared = sum(is_shared),
      n_unique = sum(is_unique),
      .groups  = "drop"
    ) %>%
    mutate(
      pct_shared = n_shared / n_total,
      pct_unique = n_unique / n_total
    )
}

loma_stats <- summarise_per_sra(lomakin_pairs)
zvy_stats  <- summarise_per_sra(zvy_pairs)

loma_partial <- loma_stats %>%
  filter(n_shared > 0, n_unique > 0) %>%
  arrange(desc(n_unique), desc(n_shared))
zvy_partial  <- zvy_stats %>%
  filter(n_shared > 0, n_unique > 0) %>%
  arrange(desc(n_unique), desc(n_shared))

keep_ids <- c("ERR6501933", "ERR6501936", "ERR6501956", "ERR6501958")

zvy_subset <- zvyagin_sra_db %>%
  filter(SRA_id %in% keep_ids)

lomakin_plus_zvy_no_meta <- bind_rows(zvy_subset, lomakin_sra_db)

fwrite(lomakin_plus_zvy_no_meta,
       "/doctorai/chiarba/add_db/lomakin/lo_plus_zvy_no_meta.tsv",
       sep = "\t", row.names = FALSE)

metadata_zvy <- fread("/doctorai/chiarba/add_db/Zvyagin/metadata_zvy.csv") %>%
  rename(SRA_id = Run)

metadata_zvy_sub <- metadata_zvy %>%
  filter(SRA_id %in% keep_ids)

desired_cols <- c(
  "SRA_id", "AGE", "Bases", "BioProject", "Broker_name", "cell_type",
  "disease", "disease_staging", "individual", "Organism_part",
  "Sample Name", "Sample_name", "sex", "SRA Study", "disease_duration"
)

metadata_keep <- metadata_zvy_sub %>%
  select(any_of(desired_cols)) %>%
  distinct(SRA_id, .keep_all = TRUE)

zvy_subset_with_meta <- zvy_subset %>%
  left_join(metadata_keep, by = "SRA_id") %>%
  mutate(study = "zvyagin")

fwrite(zvy_subset_with_meta,
       "/doctorai/chiarba/add_db/Zvyagin/zvyagin_meta_sub_uniqe.tsv",
       sep = "\t", row.names = FALSE)

metadata_lomakin <- fread("/doctorai/chiarba/add_db/lomakin/lomakin_metadata.csv") %>%
  rename(SRA_id = Run)

desired_cols_loma <- c(
  "SRA_id", "AGE", "cell_type", "common_name", "Consent", "disease",
  "Sample Name", "sex", "SRA Study", "Submitter_Id", "Library Name",
  "disease_duration", "disease_staging", "EDSS", "treatment"
)

metadata_lomakin_keep <- metadata_lomakin %>%
  select(any_of(desired_cols_loma)) %>%
  distinct(SRA_id, .keep_all = TRUE)

lomakin_with_meta <- lomakin_sra_db %>%
  left_join(metadata_lomakin_keep, by = "SRA_id")

fwrite(lomakin_with_meta,
       "/doctorai/chiarba/add_db/lomakin/lomakin_meta.tsv",
       sep = "\t", row.names = FALSE)

lomakin_plus_zvy <- bind_rows(
  zvy_subset_with_meta,
  lomakin_with_meta
) %>%
  mutate(study = coalesce(study, "lomakin"))

fwrite(lomakin_plus_zvy,
       "/doctorai/chiarba/add_db/lomakin/lomakin_plus_zvy.tsv",
       sep = "\t", row.names = FALSE)

lomakin_plus_zvy_selected <- lomakin_plus_zvy %>%
  select(
    SRA_id,
    sequence,
    junction,
    sequence_alignment,
    germline_alignment,
    cdr3,
    cdr3_aa,
    v_call,
    d_call,
    j_call,
    junction_aa,
    Organism_part,
    junction_length,
    `Sample Name`,
    treatment,
    Sample_name,
    individual,
    `Library Name`,
    study,
    cell_type
  ) %>%
  mutate(
    v_gene = sub("\\*.*", "", v_call),
    d_gene = sub("\\*.*", "", d_call),
    j_gene = sub("\\*.*", "", j_call)
  )

lomakin_plus_zvy_selected <- subset(lomakin_plus_zvy_selected, grepl("^IGH", v_gene))

lomakin_plus_zvy_selected <- lomakin_plus_zvy_selected %>%
  filter(!is.na(junction_aa) & junction_aa != "") %>%
  filter(!grepl("\\*", junction_aa)) %>%
  filter(!grepl("_",  junction_aa))

lomakin_plus_zvy_selected <- lomakin_plus_zvy_selected[grepl("^C.*W$", lomakin_plus_zvy_selected$junction_aa), ]

any(grepl("^TR", lomakin_plus_zvy_selected$v_call))
any(grepl("^TR", lomakin_plus_zvy_selected$d_gene))
any(grepl("^TR", lomakin_plus_zvy_selected$j_call))

lomakin_plus_zvy_selected <- lomakin_plus_zvy_selected %>%
  mutate(
    v_cdr3_j = paste(v_gene, junction_aa, j_gene, sep = "_")
  )

lomakin_plus_zvy_selected <- lomakin_plus_zvy_selected %>%
  mutate(
    .ms_num       = str_extract(`Library Name`, "(?<=^MS)\\d+"),
    .lo_from_lib  = if_else(!is.na(.ms_num), paste0("Lo", .ms_num), NA_character_),
    subject_number = case_when(
      individual == "MS3" ~ "Z3",
      individual == "MS5" ~ "Z5",
      TRUE                ~ .lo_from_lib
    )
  ) %>%
  select(-.ms_num, -.lo_from_lib)

lomakin_plus_zvy_selected <- lomakin_plus_zvy_selected %>%
  rename(tissue = Organism_part) %>%
  mutate(
    tissue      = "blood",
    tissue_type = case_when(tissue == "blood" ~ "extra"),
    junction_aa_length = as.numeric(junction_length) / 3,
    disease_diagnosis  = "MS"
  )

lomakin_plus_zvy_selected <- lomakin_plus_zvy_selected %>%
  filter(subject_number %in% c("Lo1", "Lo2", "Lo3", "Lo4", "Lo6", "Lo9"))

unique(lomakin_plus_zvy_selected$cell_type)

lomakin_plus_zvy_selected <- lomakin_plus_zvy_selected %>%
  filter(cell_type == "B cell")

fwrite(lomakin_plus_zvy_selected,
       "/doctorai/chiarba/add_db/lomakin/lomakin_plus_zvy_for_my_analysis.tsv",
       sep = "\t", row.names = FALSE)

#### Preprocessing HC (Ghraichy) ####

folder_path <- "/doctorai/chiarba/bcr_from_raw/HC/AIRR_HC"

tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

# Read AIRR tsv and tag with SRA_id
read_and_tag <- function(file_path) {
  sra_id <- tools::file_path_sans_ext(basename(file_path))
  df     <- fread(file_path, colClasses = "character")
  df$SRA_id <- sra_id
  df
}

HC_mixcr <- data.table::rbindlist(
  lapply(tsv_files, read_and_tag),
  use.names = TRUE,
  fill = TRUE
)

metadata_HC <- fread("/doctorai/chiarba/bcr_from_raw/HC/AIRR_HC/metadata_HC_SraRunTable.csv")

meta_subset <- metadata_HC %>%
  rename(SRA_id = Run) %>%
  select(SRA_id, sex, tissue, age, disease, `Library Name`, isolate)

HC_mixcr_meta <- HC_mixcr %>%
  left_join(meta_subset, by = "SRA_id")

fwrite(
  HC_mixcr_meta,
  "/doctorai/chiarba/bcr_from_raw/HC/AIRR_HC/HC_all_meta.tsv",
  sep = "\t",
  row.names = FALSE
)

# Select columns with standard AIRR fields used in your analysis
HC_mixcr_selected <- HC_mixcr_meta %>%
  select(
    SRA_id,
    sequence,
    junction,
    sequence_alignment,
    germline_alignment,
    cdr3,
    cdr3_aa,
    v_call,
    d_call,
    j_call,
    junction_aa,
    disease,
    tissue,
    age,
    sex,
    junction_length,
    `Library Name`,
    isolate
  ) %>%
  mutate(
    v_gene = sub("\\*.*", "", v_call),
    d_gene = sub("\\*.*", "", d_call),
    j_gene = sub("\\*.*", "", j_call)
  )

# Keep only IGH
HC_mixcr_selected <- subset(HC_mixcr_selected, grepl("^IGH", v_gene))

# Clean junction_aa
HC_mixcr_selected <- HC_mixcr_selected %>%
  filter(!is.na(junction_aa) & junction_aa != "") %>%
  filter(!grepl("\\*", junction_aa)) %>%
  filter(!grepl("_",  junction_aa))

# Keep only complete CDR3 (C...W)
HC_mixcr_selected <- HC_mixcr_selected[
  grepl("^C.*W$", HC_mixcr_selected$junction_aa),
]

# Create v_cdr3_j
HC_mixcr_selected <- HC_mixcr_selected %>%
  mutate(v_cdr3_j = paste(v_gene, junction_aa, j_gene, sep = "_"))

# Check for possible TCR contamination
any(grepl("^TR", HC_mixcr_selected$j_gene))
any(grepl("^TR", HC_mixcr_selected$d_gene))
any(grepl("^TR", HC_mixcr_selected$v_gene))

# Create subject_number
HC_mixcr_selected <- HC_mixcr_selected %>%
  rename(subject_number = isolate) %>%
  mutate(
    subject_number = str_extract(subject_number, "\\d{3}"),
    subject_number = paste0("HC_", subject_number)
  )

# Study label
HC_mixcr_selected <- HC_mixcr_selected %>%
  mutate(study = "Ghraichy")

# Fix tissue spelling
HC_mixcr_selected <- HC_mixcr_selected %>%
  mutate(tissue = dplyr::recode(tissue, "periferal blood" = "blood"))

# Rename disease column
HC_mixcr_selected <- HC_mixcr_selected %>%
  rename(disease_diagnosis = disease)

# Create junction_aa_length
HC_mixcr_selected <- HC_mixcr_selected %>%
  mutate(junction_aa_length = as.numeric(junction_length) / 3)

any(is.na(HC_mixcr_selected$junction_aa_length))

# Drop junction_length and Library Name
HC_mixcr_selected <- HC_mixcr_selected %>%
  select(-junction_length, -`Library Name`)

fwrite(
  HC_mixcr_selected,
  "/doctorai/chiarba/bcr_from_raw/HC/AIRR_HC/HC_cleaned_selectedcol_IGH.tsv",
  sep = "\t",
  row.names = FALSE
)

unique(HC_mixcr_selected$age)
length(unique(HC_mixcr_selected$age))

# Age filter (keep subjects >= 15 years)
HC_mixcr_selected_age <- HC_mixcr_selected %>%
  filter(age >= 15)

unique(HC_mixcr_selected_age$age)
length(unique(HC_mixcr_selected_age$subject_number))

fwrite(
  HC_mixcr_selected_age,
  "/doctorai/chiarba/bcr_from_raw/HC/AIRR_HC/HC_cleaned_selected_age.tsv",
  sep = "\t",
  row.names = FALSE
)

#### Build final combined MS + HC (naive) dataset ####

green  <- fread("/doctorai/chiarba/bcr_from_raw/green/AIRR/Green_cleaned_selectedcol_igh.tsv")
HC     <- fread("/doctorai/chiarba/bcr_from_raw/HC/AIRR_HC/HC_cleaned_selected_age.tsv")
pala   <- fread("/doctorai/chiarba/bcr_from_raw/pala/pala_cleaned_selectedcol_igh.tsv")

# Stern without Sanger sequences
stern  <- fread("/doctorai/chiarba/bcr_from_raw/stern/fastq_download/AIRR/stern_cleaned_selected_igh.tsv")

ramesh <- fread("/doctorai/chiarba/bcr_from_raw/ramesh_final/AIRR_format/ramesh_final_preprocessing.tsv")
lau    <- fread("/doctorai/chiarba/add_db/Laurent/laurent_naive_IgG_IgM_for_my_analysis.tsv")
lom_z  <- fread("/doctorai/chiarba/add_db/lomakin/lomakin_plus_zvy_for_my_analysis.tsv")
rus    <- fread("/doctorai/chiarba/add_db/Ruschil/ruschil_for_my_analysis.tsv")
arg    <- fread("/doctorai/chiarba/add_db/agrafiotis/agr_colse_for_my_analysis.tsv")
p_s    <- fread("/doctorai/chiarba/add_db/Perez_salvidar/IGH/perez_naive_for_my_analysis.tsv")

# Ensure SRA_id is character everywhere
green$SRA_id  <- as.character(green$SRA_id)
HC$SRA_id     <- as.character(HC$SRA_id)
stern$SRA_id  <- as.character(stern$SRA_id)
ramesh$SRA_id <- as.character(ramesh$SRA_id)
pala$SRA_id   <- as.character(pala$SRA_id)
arg$SRA_id    <- as.character(arg$SRA_id)
lau$SRA_id    <- as.character(lau$SRA_id)
p_s$SRA_id    <- as.character(p_s$SRA_id)
rus$SRA_id    <- as.character(rus$SRA_id)
lom_z$SRA_id  <- as.character(lom_z$SRA_id)

# Bind all studies (MS + HC)
MS_HC_all <- dplyr::bind_rows(
  green,
  HC,
  stern,
  ramesh,
  pala,
  arg,
  lau,
  lom_z,
  p_s,
  rus
)

unique(p_s$tretment)

# Fill missing tissue_type as extra (blood) for safety
MS_HC_all <- MS_HC_all %>%
  mutate(
    tissue_type = tidyr::replace_na(tissue_type, "extra")
  )

length(unique(MS_HC_all$subject_number[startsWith(MS_HC_all$subject_number, "HC")]))

fwrite(
  MS_HC_all,
  "/doctorai/chiarba/analysis/MS_HC_naive.tsv",
  sep = "\t",
  row.names = FALSE
)
