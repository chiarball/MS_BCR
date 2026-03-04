# 🧬 MS BCR Database 🔬


>📊 A curated B-cell receptor (BCR) repertoire database for Multiple Sclerosis research.  
>📦 The database is hosted on [Zenodo](https://zenodo.org/records/18862069) and distributed as a set of modular CSV files designed for flexible and easy-to-use downstream analysis.


---

## 📂 Repository Structure

The database is split into four CSV files, each serving a distinct purpose:

| File | Description |
|------|-------------|
| `MS_BCR_db_core.csv` | Core features for every sequence (recommended starting point) |
| `MS_BCR_db_alignment.csv` | Full alignment data, germline sequences, and per-region annotations |
| `MS_BCR_db_metadata_sparse_long.csv` | Study-specific metadata in long format (sparse fields) |
| `MS_BCR_db_miscellaneous.csv` | Additional metadata columns |

**All files share the `sequence_id` column as a primary key and can be joined as needed.**
Metadata fields with more than 50% missing values across the dataset were separated from the core table and stored in the long-format annotations file `MS_BCR_db_metadata_sparse_long.csv` to avoid sparse wide matrices. 

### Healthy Control Dataset

To facilitate direct comparison analyses, a harmonized and preprocessed healthy control dataset is also available: `MS_BCR_healthy_ONLY.csv`. The data originate from the Ghraichy et al. study and have been processed using the same pipeline applied to the MS-BCR database to ensure methodological consistency.

The dataset can be downloaded from Zenodo and corresponds to the exact control cohort used in the analyses presented in our manuscript.

---

## File Descriptions

### `MS_BCR_db_core.csv`
The main file for most analyses. Contains one row per sequence with the following fields:

| Column | Description |
|--------|-------------|
| `sequence_id` | Unique sequence identifier (primary key) |
| `clone_id` | Clonal cluster assignment |
| `study` | Source study |
| `subject_number` | Donor identifier |
| `disease_diagnosis` | Clinical diagnosis (e.g. `Healthy`) |
| `sex` | Donor sex |
| `tissue` / `tissue_type` | Tissue of origin |
| `locus` | Receptor locus (`IGH`, `IGK`, `IGL`) |
| `cell_subtype` | B-cell phenotype (e.g. `Naive B cell`, `Memory B cell`, `Plasmablast`) |
| `v_gene`, `d_gene`, `j_gene` | VDJ gene calls (allele stripped) |
| `junction_aa` | Junction amino acid sequence |
| `junction_length_aa` | Junction length in amino acids |
| `productive` | Whether the sequence is productive |

---

### `MS_BCR_db_alignment.csv`
Contains raw sequences, alignment coordinates, CIGAR strings, CDR/FWR region sequences, and alignment scores. Useful for structural or mutation analyses.

---

### `MS_BCR_db_miscellaneous.csv`
Contains time_point,treatment,BioProject,BioSample,Sample Name,EDSS,age,Bases,Bytes. Miscellaneous metadata.

---

### `MS_BCR_db_metadata_sparse_long.csv`
Study-specific metadata (e.g. time points, treatment, disease scores) stored in a **long format**. Schema:

| Column | Description |
|--------|-------------|
| `sequence_id` | Primary key |
| `field` | Metadata field name |
| `value` | Metadata value |

---

## 🚀 Quick Start

### Load the core database
```r
library(data.table)

core <- fread("MS_BCR_db_core.csv")
```

### Join with alignment data
```r
alignment <- fread("MS_BCR_db_alignment.csv")

db_full <- merge(core, alignment, by = "sequence_id", all.x = TRUE)
```

### Merge study-specific metadata
```r
library(tidyr)
library(dplyr)

annotations <- fread("MS_BCR_db_metadata_sparse_long.csv")

# Convert to wide format
annotations_wide <- annotations %>%
  pivot_wider(names_from = field, values_from = value)

# Join into core
core_annotated <- core %>%
  left_join(annotations_wide, by = c("sequence_id", "study", "subject_number"))
```

---


## Notes and database preprocessing details

- All gene calls have allele information stripped (e.g. `IGHV1-69*01` → `IGHV1-69`).
- The metadata file uses long format intentionally — most study-specific fields are highly sparse and this format avoids large empty dataframes.
- `sequence_id` is globally unique across studies and can be safely used as a join key across all files.

Data Processing
Cell subtype labels were standardized across studies to harmonize inconsistent nomenclature (e.g. "Memory B lymphocyte" → "Memory B cell").  The receptor locus (IGH, IGK, IGL) was inferred from the V gene call. Metadata fields with more than 50% missing values across the dataset were separated from the core table and stored in a long-format annotations file to avoid sparse wide matrices. 


---

## Citation

If you use this database, please cite:

> [Citation placeholder]

---

## Contact

For questions or issues, please open a GitHub issue or contact [your email/contact].