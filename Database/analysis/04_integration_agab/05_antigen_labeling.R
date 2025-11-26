#!/usr/bin/env Rscript

# Merge Fisher antigen-level results (lev0/1/2) and keep only significant rows

  library(readr)
  library(dplyr)

# --- INPUT BASE DIRECTORY (clean version) ---
base_dir <- "/doctorai/chiarba/AbAg_database/clean"

# --- INPUT FILES (from step 3 / Python Fisher) ---
files <- c(
  lev0 = file.path(base_dir, "fisher_results_annot1_antigen_lev0.csv"),
  lev1 = file.path(base_dir, "fisher_results_annot1_antigen_lev1.csv"),
  lev2 = file.path(base_dir, "fisher_results_annot1_antigen_lev2.csv")
)

# --- READ EACH LEVEL AND ADD 'lev' COLUMN ---
dfs <- lapply(names(files), function(nm) {
  df <- read_csv(files[[nm]], show_col_types = FALSE, guess_max = 100000)
  df$lev <- nm
  df
})

# Row-bind all levels
all_levels <- bind_rows(dfs)

# Sanity check for p_adj
if (!"p_adj" %in% names(all_levels)) {
  stop("Column 'p_adj' not found. Please check your CSV headers.")
}

# Ensure p_adj is numeric
all_levels <- all_levels %>%
  mutate(p_adj = as.numeric(p_adj))

# --- FILTER SIGNIFICANT RESULTS AND SORT ---
sig <- all_levels %>%
  filter(!is.na(p_adj) & p_adj < 0.05) %>%
  arrange(lev, p_adj)

# --- OUTPUT (clean, non-overwriting name) ---
out_path <- file.path(base_dir, "fisher_results_annot1_antigen_all_levels_sig.csv")
write_csv(sig, out_path)

# --- SUMMARY ON STDOUT ---
cat(
  "Merged rows:", nrow(all_levels), "\n",
  "Significant rows (p_adj < 0.05):", nrow(sig), "\n",
  "Saved to:", out_path, "\n"
)

```

```{r step5 annotation taxonomy group}

#!/usr/bin/env Rscript

# Step 5: classify antigen groups into broad origin classes (label)
# - Input:  fisher_results_annot1_antigen_all_levels_sig.csv
# - Output: fisher_labeled_interestMS1.tsv
# Code comments in English only.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
})

# ---------------------- I/O PATHS ----------------------

base_dir  <- "/doctorai/chiarba/AbAg_database/clean"
in_path   <- file.path(base_dir, "fisher_results_annot1_antigen_all_levels_sig.csv")
out_path  <- file.path(base_dir, "fisher_labeled_interestMS1.tsv")

# ---------------------- LOAD DATA ----------------------

sig_ms <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)

if (!"group" %in% names(sig_ms)) {
  stop("Column 'group' not found in input file. Check Fisher output.")
}

# ---------------------- BUILD GROUP TABLE ----------------------

# Unique group names to classify
gvec <- unique(as.character(sig_ms$group))
df   <- tibble(group = gvec, group_lc = tolower(gvec))

# ---------------------- REGEX RULES (MERGED VERSION) ----------------------

# 1) Toxins
rx_toxin <- c(
  "\\btoxin\\b", "neurotoxin", "exotoxin", "endotoxin",
  "alpha-?hemolysin", "aureolysin", "dermonecrotic", "vaginolysin", "botulinum"
)

# 2) Parasites (protozoa/helminths)
rx_parasite <- c(
  "plasmodium|merozoite|ookinete|circumsporozoite",
  "toxoplasma|\\bsag1\\b", "leishmania", "trypanosom|\\bvsg\\b",
  "schisto", "echinococcus", "taenia", "onchocerca", "giardia", "entamoeba", "trichomon"
)

# 3) Viral (more inclusive merged list)
rx_viral <- c(
  # Explicit Spike first
  "\\bspike glycoprotein\\b",
  # EBV / CMV
  "epstein[- ]?barr|\\bebv\\b",
  "cytomegalovirus|\\bhcmv\\b|\\bcmv\\b",
  # Other viruses
  "adenovirus|herpes(virus)?",
  "influenza|hemagglutinin|neuraminidase",
  "\\bhiv\\b|envelope glycoprotein gp\\d+|gp160",
  "papillomavirus|\\bhpv\\b",
  "hepatitis [abcdeg]|\\bhbv\\b|\\bhcv\\b",
  "coronavirus|sars[-_ ]?cov\\w*",
  # Generic viral ORFs and structural proteins
  "\\borf\\d+\\b|orf3a|orf7a|orf9b",
  "\\bnucleoprotein\\b|\\bcapsid protein\\b|\\benvelope protein\\b|\\bmatrix protein\\b",
  "\\breplicase\\b|polyprotein",
  # AAV family
  "\\baav\\b|adeno[- ]associated virus",
  "assembly activating protein|\\baap\\b",
  "protein rep\\s?(40|52|68|78)\\b|\\brep(40|52|68|78)\\b|\\bvp1\\b",
    "immediate early protein ie1|immediate early protein",
  "\\bul\\d+\\b",
  "non-structural protein"
)

# 4) Bacterial (merged list, includes Akkermansia)
rx_bacterial <- c(
  "akkermansia", "staphyl", "streptoc", "pseudomonas", "mycobacter", "clostrid",
  "escherichia|\\becoli\\b", "salmonella", "listeria", "neisseria", "borrelia",
  "brucella", "vibrio", "campylobacter", "bacteroides",
  "outer membrane", "lipopolysaccharide|\\blps\\b"
)

# 5) Plant
rx_plant <- c(
  "pollen|allergen", "\\bphl p\\b", "bet v 1", "\\blol p\\b",
  "arabidopsis", "\\bplant\\b"
)

# 6) Human (core) – broad but with viral exclusions
rx_human_core <- c(
  "\\bhla\\b|histocompatibility",
  "\\bcd\\d+\\b",
  "interleukin|\\bil-?\\d",
  "tumor necrosis factor|\\btnf\\b",
  "chemokine|\\b(?:cxcl|ccl|xcl|cx3cl)\\d+\\b",
  "immunoglobulin|heavy constant|kappa\\b|lambda\\b|fc receptor|fcrn|mhc",
  "receptor|integrin|kinase|phosphatase|transporter|enzyme|growth factor|hormone|collagen|matrix|protease|glycoprotein",
  "neurotensin receptor|adrenergic|dopamine|serotonin|acetylcholine receptor",
  "albumin|ferritin|transferrin|actinin|myosin|dystroph|sarcoglycan|utrophin|laminin",
  "beta-?amyloid|\\btau\\b|alpha-?synuclein|prion|psen\\d|pink1",
  "insulin|leptin|ghrelin|calcitonin|glucagon|thyroid|prolactin",
  "\\bjak\\b|\\bstat\\b|mapk|erk|\\bakt\\b|raf\\b|pi3k|mtor|tgf-?beta|bmp\\b|smad\\b|\\bwnt\\b|frizzled|\\bnotch\\b",
  "carboxylesterase 1\\b|\\blysozyme c\\b"
)

# Words that, if present, should block a "human" label
viral_exclusion <- paste(
  c(
    "spike glycoprotein",
    "capsid",
    "nucleoprotein",
    "\\borf\\d+\\b",
    "polyprotein",
    "replicase",
    "aav",
    "gp\\d+",
    "immediate early protein",
    "\\bul\\d+\\b",
    "non-structural protein"
  ),
  collapse = "|"
)


# ---------------------- CLASSIFIER ----------------------

classify_group <- function(x) {
  # x is lowercased group name
  pick_any <- function(rx_vec) {
    any(str_detect(x, regex(paste(rx_vec, collapse = "|"), ignore_case = TRUE)))
  }
  if (pick_any(rx_toxin))     return("toxin")
  if (pick_any(rx_parasite))  return("parasite")
  if (pick_any(rx_viral))     return("viral")
  if (pick_any(rx_bacterial)) return("bacterial")
  # human only if no clear viral markers
  if (!str_detect(x, regex(viral_exclusion, ignore_case = TRUE)) &&
      pick_any(rx_human_core)) return("human")
  if (pick_any(rx_plant))     return("plant")
  return("other")
}

df <- df %>%
  mutate(broad_origin = vapply(group_lc, classify_group, character(1)))

# ---------------------- SIMPLE OVERRIDE MAP (STRING MATCH) ----------------------

override_map <- c(
  # explicit fixes
  "lysozyme c"         = "human",
  "spike glycoprotein" = "viral",
  "protein rep52"      = "viral",
  # EBV/CMV/Akkermansia (reinforce classification)
  "epstein-barr"       = "viral",
  "ebv"                = "viral",
  "cytomegalovirus"    = "viral",
  "akkermansia"        = "bacterial",
  # human enzymes that may look generic
  "liver carboxylesterase 1" = "human"
)

for (k in names(override_map)) {
  df$broad_origin[str_detect(df$group_lc, fixed(k, ignore_case = TRUE))] <- override_map[[k]]
}

# ---------------------- QUICK SUMMARY (BEFORE MANUAL ANTIGEN OVERRIDES) ---------

cat("Unique groups per class (before manual antigen overrides):\n")
print(table(df$broad_origin))

# ======================================================================
# MANUAL ANTIGEN-LEVEL OVERRIDES (based on inspection of 'other' groups)
# ======================================================================

# --- Block 1: main human / viral list --------------------------------

human_list <- c(
  "Programmed cell death 1 ligand 1",
  "Angiotensinogen",
  "Beta-synuclein",
  "Cytochrome c, somatic",
  "GTPase KRas",
  "Proprotein convertase subtilisin/kexin type 9",
  "Epididymal-specific lipocalin-6",
  "Dehydrogenase/reductase SDR family member 11",
  "Syntaxin-6",
  "Syntaxin-5",
  "Sterol regulatory element-binding protein cleavage-activating protein",
  "Sodium-dependent phosphate transport protein 4",
  "Neutral ceramidase",
  "Solute carrier family 22 member 6",
  "Ceramide synthase 1",
  "Intersectin-1",
  "Syntaxin-17",
  "Solute carrier family 22 member 11",
  "Cathepsin B-like cysteine proteinase",
  "Lymphocyte antigen 6S",
  "Rho GTPase-activating protein 35",
  "Vesicle-associated membrane protein 2",
  "Rhombotin-2",
  "Dipeptidyl peptidase 3",
  "Gamma-synuclein",
  "Calreticulin-3",
  "Peptidyl-prolyl cis-trans isomerase FKBP4",
  "Centrosomal protein of 57 kDa",
  "Aldo-keto reductase family 1 member D1",
  "Proteasome subunit beta type-4",
  "MARCKS-related protein",
  "Protein ADM2",
  "Cytochrome P450 11B2, mitochondrial",
  "Ceramide synthase 5",
  "Protein disulfide-isomerase A3",
  "Mucosal addressin cell adhesion molecule 1",
  "GTPase HRas",
  "Aflatoxin B1 aldehyde reductase member 2",
  "Aldo-keto reductase family 1 member B1",
  "Aldo-keto reductase family 1 member A1",
  "Aldo-keto reductase family 7 member A3",
  "Aflatoxin B1 aldehyde reductase member 4",
  "Cadherin-6",
  "Peptidyl-Asp metalloendopeptidase",
  "Glutathione S-transferase A3",
  "Cytochrome P450 2W1",
  "Pro-neuregulin-2, membrane-bound isoform",
  "Protein Tob2",
  "Mesothelin",
  "Synphilin-1",
  "Bcl-2 homologous antagonist/killer",
  "Pro-adrenomedullin",
  "Glomulin",
  "Cryptic family protein 1B",
  "Cryptic protein",
  "Glucocorticoid modulatory element-binding protein 1",
  "Huntingtin",
  "Cytochrome P450 2A6",
  "Nicotinamide phosphoribosyltransferase",
  "Programmed cell death protein 1",
  "Rod cGMP-specific 3',5'-cyclic phosphodiesterase subunit alpha",
  "Rod cGMP-specific 3',5'-cyclic phosphodiesterase subunit beta",
  "High affinity cGMP-specific 3',5'-cyclic phosphodiesterase 9A",
  "NKG2-D type II integral membrane protein",
  "Carbonic anhydrase 2",
  "Cytotoxic T-lymphocyte protein 4",
  "Cathepsin B",
  "Glypican-3",
  "Programmed cell death 1 ligand 2",
  "Poly [ADP-ribose] polymerase tankyrase-2",
  "Poly [ADP-ribose] polymerase tankyrase-1",
  "Claudin-18",
  "Delta-like protein 3",
  "V-set domain-containing T-cell activation inhibitor 1",
  "Claudin-6",
  "Melanoma-associated antigen B18",
  "Melanoma-associated antigen 12",
  "Melanoma-associated antigen 1",
  "Melanoma-associated antigen 4",
  "Melanoma-associated antigen C2",
  "Protein LSM14 homolog A",
  "Coagulation factor VIII",
  "Neutrophil elastase",
  "Melanoma-associated antigen D4",
  "Melanoma-associated antigen B4",
  "Melanoma-associated antigen C3",
  "Melanoma-associated antigen 6",
  "Coagulation factor X",
  "Melanoma-associated antigen 3",
  "Cell adhesion molecule CEACAM5",
  "Lymphocyte antigen 75",
  "Caspase recruitment domain-containing protein 14",
  "Caspase recruitment domain-containing protein 16",
  "Caspase recruitment domain-containing protein 18",
  "Caspase recruitment domain-containing protein 19",
  "Putative caspase recruitment domain-containing protein 17P"
)

viral_list <- c(
  "Envelope small membrane protein",
  "Membrane protein",
  "RNA-directed RNA polymerase L"
)

df$broad_origin[df$group %in% human_list] <- "human"
df$broad_origin[df$group %in% viral_list] <- "viral"

# --- Block 2: reclassification for the next chunk -------------------

human_list2 <- c(
  "Putative caspase recruitment domain-containing protein 17P",
  "Caspase recruitment domain-containing protein 9",
  "Caspase recruitment domain-containing protein 11",
  "Recoverin",
  "Caspase recruitment domain-containing protein 10",
  "Coagulation factor IX",
  "Mucin-1",
  "Protein CD300H",
  "Alanine--tRNA ligase, cytoplasmic",
  "Glutamate carboxypeptidase 2",
  "Protein S100-A9",
  "Stimulator of interferon genes protein",
  "Mucin-16",
  "TOX high mobility group box family member 2",
  "Basigin",
  "Somatotropin",
  "Prolyl endopeptidase FAP",
  "Sialic acid-binding Ig-like lectin 15",
  "Breakpoint cluster region protein",
  "C-type lectin domain family 12 member A",
  "Induced myeloid leukemia cell differentiation protein Mcl-1",
  "Nectin-4",
  "HAUS augmin-like complex subunit 3",
  "C2 domain-containing protein 3",
  "Lysine-specific histone demethylase 1A",
  "Apoptosis regulator Bcl-2",
  "Annexin A1",
  "Kit ligand",
  "Lymphocyte antigen 6E",
  "Progranulin",
  "Metastasis-suppressor KiSS-1",
  "GRIP and coiled-coil domain-containing protein 2",
  "Granulocyte-macrophage colony-stimulating factor",
  "Ankyrin repeat domain-containing protein 36B",
  "Neural cell adhesion molecule L1",
  "Zinc finger and BTB domain-containing protein 32",
  "SLAM family member 7",
  "Tenascin-N",
  "Basal cell adhesion molecule",
  "Sortilin",
  "Fibrinogen alpha chain",
  "Interferon gamma",
  "Tenascin",
  "Fibrinogen beta chain",
  "Elongation factor 1-delta",
  "Inducible T-cell costimulator",
  "Fibrinogen gamma chain",
  "TAR DNA-binding protein 43",
  "Cell adhesion molecule CEACAM8",
  "C-type lectin domain family 12 member B",
  "14-3-3 protein theta",
  "Lymphocyte activation gene 3 protein",
  "GRIP and coiled-coil domain-containing protein 1",
  "Obscurin",
  "Obscurin-like protein 1",
  "Telethonin",
  "F-box only protein 9",
  "Wilms tumor protein",
  "PHD finger protein 20",
  "Melanoma-associated antigen 2",
  "Splicing factor Cactin",
  "Mucin-17",
  "NudC domain-containing protein 1",
  "Melanoma-associated antigen D1",
  "Glypican-2",
  "Adenosine deaminase 2",
  "Cell adhesion molecule CEACAM6",
  "Titin",
  "Signal transducer and activator of transcription 3",
  "Tumor-associated calcium signal transducer 2",
  "Histone-lysine N-methyltransferase EHMT1",
  "Neural cell adhesion molecule 1",
  "Islet amyloid polypeptide",
  "Zinc finger CCCH domain-containing protein 14",
  "Apoptosis regulatory protein Siva",
  "C-type lectin domain family 7 member A",
  "Stomatin-like protein 2, mitochondrial",
  "Coagulation factor XI",
  "Endosialin",
  "Heat shock 70 kDa protein 14",
  "Delta-like protein 1",
  "Platelet endothelial cell adhesion molecule",
  "ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 2",
  "Interferon beta",
  "Cell adhesion molecule CEACAM21",
  "SH2 domain-containing protein 2A",
  "Transmembrane and coiled-coil domain protein 3",
  "E3 ubiquitin-protein ligase XIAP",
  "Thrombospondin-2",
  "NEDD4 family-interacting protein 1",
  "Roundabout homolog 2",
  "p53 apoptosis effector related to PMP-22",
  "B-cell linker protein",
  "Talin-2",
  "VA11.1 protein",
  "E3 ubiquitin-protein ligase NHLRC1"
)

toxin_list2 <- c(
  "Alpha-cobratoxin"
)

viral_list2 <- c(
  "Non-structural protein ORF4a",
  "Non-structural protein ORF4b"
)

bacterial_list2 <- c(
  "Large ribosomal subunit protein uL22",
  "Integration host factor"
)

df$broad_origin[df$group %in% human_list2]     <- "human"
df$broad_origin[df$group %in% toxin_list2]     <- "toxin"
df$broad_origin[df$group %in% viral_list2]     <- "viral"
df$broad_origin[df$group %in% bacterial_list2] <- "bacterial"

# --- Block 3 --------------------------------------------------------

human_list3 <- c(
  "E3 ubiquitin-protein ligase NHLRC1",
  "E3 ubiquitin-protein ligase NEDD4",
  "Histone H2B type 1",
  "Alpha-amylase 2B",
  "DNA-directed RNA polymerase II subunit RPB1",
  "Thioredoxin",
  "Metallothionein",
  "F-box DNA helicase 1",
  "Protein S100-A3",
  "Major vault protein",
  "Myotubularin",
  "Sphingomyelinase",
  "V(D)J recombination-activating protein 2",
  "V(D)J recombination-activating protein 1",
  "Protein HEG homolog 1",
  "PALM2-AKAP2 fusion protein",
  "Rho GTPase-activating protein 7",
  "Bromodomain-containing protein 4",
  "15-hydroxyprostaglandin dehydrogenase [NAD(+)]",
  "Peroxynitrite isomerase 1",
  "Hematopoietic prostaglandin D synthase",
  "Laforin",
  "Plexin domain-containing protein 1",
  "Plexin domain-containing protein 2",
  "Aldo-keto reductase family 1 member C3",
  "Muscleblind-like protein 1",
  "Interferon alpha-1/13",
  "Sialic acid-binding Ig-like lectin 5",
  "Interferon alpha-8",
  "Interferon alpha-7",
  "Cell adhesion molecule CEACAM4",
  "Interferon alpha-17",
  "Interferon alpha-14",
  "Sialic acid-binding Ig-like lectin 7",
  "Interferon alpha-21",
  "Claudin-1",
  "Chondroitin sulfate proteoglycan 4",
  "Desmoglein-3",
  "Sialic acid-binding Ig-like lectin 9",
  "Interferon alpha-4",
  "Cell surface A33 antigen",
  "B- and T-lymphocyte attenuator",
  "Sialic acid-binding Ig-like lectin 6",
  "Cell adhesion molecule CEACAM19",
  "Cell adhesion molecule CEACAM3",
  "Cell adhesion molecule CEACAM20",
  "Cell adhesion molecule CEACAM16",
  "Cell adhesion molecule CEACAM7",
  "Agrin",
  "Differentially expressed in FDCP 6 homolog",
  "Neuroligin-3",
  "Lymphocyte antigen 6 complex locus protein G6d",
  "Inactive rhomboid protein 2",
  "Regulator of MON1-CCZ1 complex",
  "Heparanase",
  "Sprouty-related, EVH1 domain-containing protein 1",
  "Inactive heparanase-2",
  "Anterior gradient protein 2 homolog",
  "MICOS complex subunit MIC19",
  "Glycosylphosphatidylinositol-anchored high density lipoprotein-binding protein",
  "BRCA2 and CDKN1A-interacting protein",
  "MICOS complex subunit MIC60",
  "Prostaglandin-H2 D-isomerase",
  "G1/S-specific cyclin-D3",
  "Cell adhesion molecule CEACAM1",
  "Retinoic acid early transcript 1E",
  "UL-16 binding protein 5",
  "UL16-binding protein 6",
  "UL16-binding protein 2",
  "UL16-binding protein 3",
  "UL16-binding protein 1",
  "Serine--tRNA ligase, cytoplasmic",
  "TP53-binding protein 1",
  "Carbonyl reductase [NADPH] 1",
  "Hemojuvelin",
  "Sestrin-1",
  "C2 calcium-dependent domain-containing protein 4A",
  "Plexin-A1",
  "Complement decay-accelerating factor",
  "Retinoblastoma-associated protein",
  "Protein FAM83B",
  "Ankyrin repeat domain-containing protein 13A",
  "WW domain binding protein VOPP1",
  "Quinone oxidoreductase PIG3",
  "Probetacellulin",
  "Ankyrin repeat domain-containing protein 13B",
  "Ras-related protein Rab-11A"
)

bacterial_list3 <- c(
  "DNA-directed RNA polymerase subunit beta'",
  "Type III secretion system protein",
  "Multicopper oxidase MmcO",
  "Probable arabinosyltransferase B",
  "Uncharacterized protein MG281"
)

plant_list3 <- c(
  "Cytokinin riboside 5'-monophosphate phosphoribohydrolase"
)

viral_list3 <- c(
  "Early E3B 10.4 kDa protein",
  "Nuclear egress protein 2",
  "F protein"
)

parasite_list3 <- c(
  "Micronemal protein 1",
  "Micronemal protein 3",
  "Micronemal protein 4",
  "Micronemal protein 6",
  "Major surface antigen p30"
)

df$broad_origin[df$group %in% human_list3]      <- "human"
df$broad_origin[df$group %in% bacterial_list3]  <- "bacterial"
df$broad_origin[df$group %in% plant_list3]      <- "plant"
df$broad_origin[df$group %in% viral_list3]      <- "viral"
df$broad_origin[df$group %in% parasite_list3]   <- "parasite"

# --- Block 4 --------------------------------------------------------

human_list4 <- c(
  "Ras-related protein Rab-11A",
  "Epigen",
  "Anoctamin-9",
  "Tumor protein p53-inducible nuclear protein 1",
  "Tumor protein p73",
  "Epithelial membrane protein 2",
  "Tissue factor",
  "Epithelial cell adhesion molecule",
  "Pancreatic alpha-amylase",
  "Complement factor D",
  "Neuronal regeneration-related protein",
  "Galectin-3-binding protein",
  "Lysosomal alpha-glucosidase",
  "Gap junction beta-5 protein",
  "Protein unc-119 homolog A",
  "Telomerase reverse transcriptase",
  "Anoctamin-3",
  "Melanoma-associated antigen 11",
  "Fatty acid-binding protein, adipocyte",
  "Heat shock factor-binding protein 1",
  "Protein Cripto",
  "Homeobox protein DLL-3",
  "Ig-like domain-containing protein",
  "Maltase-glucoamylase",
  "Apoptosis-stimulating of p53 protein 2",
  "Lipoprotein lipase",
  "Cholinesterase",
  "Plasminogen activator inhibitor 2",
  "Plasma kallikrein",
  "Repulsive guidance molecule B",
  "Neogenin",
  "Prostate-specific antigen",
  "Granzyme M",
  "Prelamin-A/C",
  "Gastrin",
  "Roundabout homolog 4",
  "Amphiregulin",
  "Protein Mdm4",
  "Reticulon-4",
  "Syntaxin-4",
  "Repulsive guidance molecule A",
  "Zinc fingers and homeoboxes protein 2",
  "ProSAAS",
  "Signal transducer and activator of transcription 1-alpha/beta",
  "Tripartite motif-containing protein 26",
  "Lysine-specific demethylase 5A",
  "Set1/Ash2 histone methyltransferase complex subunit ASH2",
  "Tryptase gamma",
  "Tryptase alpha/beta-1",
  "G1/S-specific cyclin-D1",
  "Tryptase beta-2",
  "Protein SET",
  "Histone H2A.Z",
  "Gastric inhibitory polypeptide",
  "Histone PARylation factor 1",
  "Acidic leucine-rich nuclear phosphoprotein 32 family member B",
  "Tissue-type plasminogen activator",
  "Acidic leucine-rich nuclear phosphoprotein 32 family member A",
  "Oncostatin-M",
  "Complement component C6",
  "Endoplasmic reticulum chaperone BiP",
  "Bone morphogenetic protein 4",
  "Histone RNA hairpin-binding protein",
  "Histone acetyltransferase type B catalytic subunit",
  "Lysine-specific histone demethylase 2",
  "Histone chaperone ASF1A",
  "Histone-binding protein RBBP7",
  "L-lactate dehydrogenase B chain",
  "Bone morphogenetic protein 7",
  "Leukosialin",
  "Tryptophan--tRNA ligase, cytoplasmic",
  "Mucin-5AC",
  "Coagulation factor V",
  "Mitochondrial import inner membrane translocase subunit TIM16",
  "Recombining binding protein suppressor of hairless",
  "Harmonin",
  "Prostaglandin E synthase 2",
  "Bone morphogenetic protein 2",
  "Protein OPG161",
  "Tryptase delta",
  "LIM and senescent cell antigen-like-containing domain protein 1",
  "CUB domain-containing protein 1",
  "ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 1",
  "Delta-like protein 4",
  "Gasdermin-C",
  "1-phosphatidylinositol 4,5-bisphosphate phosphodiesterase gamma-1",
  "Caspase-2",
  "Tomoregulin-2",
  "Tenascin-X",
  "Tenascin-R",
  "E3 ubiquitin-protein ligase RNF123",
  "L-selectin",
  "Sialomucin core protein 24",
  "Signal-transducing adaptor protein 1",
  "THO complex subunit 5"
)

viral_list4 <- c(
  "Protein BNLF2a",
  "Accessory protein p12I",
  "Protein Tax-1",
  "Immediate early protein IE1",
  "Structural protein",
  "ICP47 protein"
)

df$broad_origin[df$group %in% human_list4] <- "human"
df$broad_origin[df$group %in% viral_list4] <- "viral"

# --- Block 5 --------------------------------------------------------

human_list5 <- c(
  "THO complex subunit 5",
  "Signal-transducing adaptor protein 2",
  "Flotillin-2",
  "Growth/differentiation factor 2",
  "Inhibin beta A chain",
  "Inhibin beta B chain",
  "E3 UFM1-protein ligase 1",
  "Lymphotoxin-beta",
  "Beta-klotho",
  "Beta-2-microglobulin",
  "CKLF-like MARVEL transmembrane domain-containing protein 4",
  "Carbonic anhydrase 9",
  "CKLF-like MARVEL transmembrane domain-containing protein 6",
  "Calreticulin",
  "Krueppel-like factor 15",
  "CCN family member 1",
  "C-type lectin domain family 14 member A",
  "Small integral membrane protein 20",
  "5'-nucleotidase",
  "Transmembrane protein 176B",
  "Cell division cycle protein 20 homolog B",
  "Lymphotoxin-alpha",
  "Macrophage colony-stimulating factor 1",
  "Growth/differentiation factor 8",
  "Follistatin",
  "Transmembrane protein 176A",
  "Lactase-like protein",
  "Proepiregulin",
  "C-type lectin domain family 4 member C",
  "Poly [ADP-ribose] polymerase 2",
  "Myc proto-oncogene protein",
  "Klotho",
  "Membrane-spanning 4-domains subfamily A member 7",
  "Membrane-spanning 4-domains subfamily A member 3",
  "Adenomatous polyposis coli protein",
  "Poly [ADP-ribose] polymerase 1",
  "Target of Myb1 membrane trafficking protein",
  "Membrane-spanning 4-domains subfamily A member 4A",
  "Membrane-spanning 4-domains subfamily A member 6A",
  "F-box only protein 38",
  "Xaa-Pro dipeptidase",
  "Butyrophilin subfamily 2 member A1",
  "Embryonic growth/differentiation factor 1",
  "High mobility group protein B1",
  "Tripartite motif-containing protein 64",
  "Beta-1-syntrophin",
  "Mediator of RNA polymerase II transcription subunit 1",
  "Myocyte-specific enhancer factor 2D",
  "Histone acetyltransferase p300",
  "Mediator of RNA polymerase II transcription subunit 31",
  "Mediator of RNA polymerase II transcription subunit 10",
  "Mediator of RNA polymerase II transcription subunit 13",
  "Glycogen [starch] synthase, muscle",
  "Alpha-1-syntrophin",
  "Mediator of RNA polymerase II transcription subunit 14",
  "Dystrobrevin alpha",
  "Mediator of RNA polymerase II transcription subunit 18",
  "Mediator of RNA polymerase II transcription subunit 19",
  "Mediator of RNA polymerase II transcription subunit 26",
  "Dystroglycan 1",
  "Mediator of RNA polymerase II transcription subunit 20",
  "Pro-neuregulin-1, membrane-bound isoform",
  "Growth/differentiation factor 7",
  "Bone morphogenetic protein 3",
  "Growth/differentiation factor 10",
  "Bone morphogenetic protein 8A",
  "Double homeobox protein 4C",
  "Protein mono-ADP-ribosyltransferase TIPARP",
  "Poly(ADP-ribose) glycohydrolase",
  "Desmin",
  "Poly [ADP-ribose] polymerase tankyrase",
  "Eosinophil peroxidase",
  "Double homeobox protein 4",
  "Poly [ADP-ribose] polymerase",
  "Mediator of RNA polymerase II transcription subunit 23",
  "Left-right determination factor 1",
  "Prostate stem cell antigen",
  "Growth/differentiation factor 11",
  "Cadherin-1",
  "Pentraxin-related protein PTX3",
  "Frataxin, mitochondrial",
  "NKG2-A/NKG2-B type II integral membrane protein",
  "Butyrophilin subfamily 3 member A1",
  "Bone morphogenetic protein 5",
  "Erythropoietin",
  "Eukaryotic translation initiation factor 4E",
  "Polycomb protein EED",
  "Mannosyl-oligosaccharide 1,2-alpha-mannosidase IA",
  "EEF1A lysine methyltransferase 3",
  "E3 ubiquitin-protein ligase UBR5",
  "Sperm-associated microtubule inner protein 4",
  "Cellular tumor antigen p53",
  "NACHT, LRR and PYD domains-containing protein 1",
  "Histone-lysine N-methyltransferase 2A",
  "Left-right determination factor 2",
  "TOM1-like protein 2",
  "Polycomb protein SUZ12",
  "Target of EGR1 protein 1",
  "Interferon omega-1",
  "Protein Tob1",
  "Forkhead box protein H1"
)

df$broad_origin[df$group %in% human_list5] <- "human"

# --- Block 6 --------------------------------------------------------

human_list6 <- c(
  "Forkhead box protein H1",
  "Zinc finger protein Eos",
  "Zinc finger protein Helios",
  "Protein ITPRID2",
  "C-type lectin domain family 9 member A",
  "Hsp90 co-chaperone Cdc37",
  "Hsp90 co-chaperone Cdc37-like 1",
  "C-type lectin domain family 6 member A",
  "C-type lectin domain family 11 member A",
  "Follistatin-related protein 3",
  "C-type lectin domain family 4 member E",
  "C-type lectin domain family 2 member B",
  "Ectonucleoside triphosphate diphosphohydrolase 2",
  "Complement C3",
  "Cytochrome P450 2C18",
  "Galactoside alpha-(1,2)-fucosyltransferase 1",
  "CCN family member 2",
  "C-type lectin domain family 4 member A",
  "Cytochrome P450 2D6",
  "C-type lectin domain family 4 member D",
  "Cytochrome P450 2C19",
  "Dexamethasone-induced protein",
  "Cytochrome P450 3A43",
  "Cytochrome P450 1A2",
  "Dexamethasone-induced Ras-related protein 1",
  "Proteasome subunit beta type-1",
  "Proteasome subunit beta type-5",
  "F-box only protein 32",
  "Chitinase domain-containing protein 1",
  "Claudin-9",
  "Dickkopf-related protein 1",
  "Cytochrome P450 2B6",
  "Cytochrome P450 2C9",
  "Histone-lysine N-methyltransferase EZH2",
  "Cytochrome P450 3A5",
  "Cytochrome P450 2C8",
  "Heat shock protein 75 kDa, mitochondrial",
  "Hyaluronidase-2",
  "Cell surface hyaluronidase CEMIP2",
  "Hyaluronidase-3",
  "Hyaluronidase-1",
  "Hyaluronidase-4",
  "Glycophorin-B",
  "Kinesin-like protein KIF11",
  "SH3 domain-binding protein 4",
  "Protein 4.2",
  "Band 3 anion transport protein",
  "Transmembrane protein PVRIG",
  "Protein S100-A5",
  "AN1-type zinc finger protein 6",
  "Somatostatin",
  "Neutrophil gelatinase-associated lipocalin",
  "Hyaluronidase PH-20",
  "Amyloid beta precursor protein binding family B member 2",
  "Cytochrome P450 3A4",
  "Caspase-7",
  "TM2 domain-containing protein 1",
  "Caspase-1",
  "Properdin",
  "ATP-dependent translocase ABCB1",
  "Phosphoprotein associated with glycosphingolipid-enriched microdomains 1",
  "Baculoviral IAP repeat-containing protein 5",
  "DNA-directed RNA polymerase III subunit RPC9",
  "Neprilysin",
  "Mucin-like protein 1",
  "C-type lectin domain family 4 member K",
  "C-type lectin domain family 4 member M",
  "E3 ubiquitin-protein ligase Mdm2",
  "Glycophorin-A",
  "Nectin-2",
  "Membrane cofactor protein",
  "Apolipoprotein B-100",
  "Signal-regulatory protein gamma",
  "Amyloid-beta precursor protein",
  "Syndecan-1",
  "Galactoside alpha-(1,2)-fucosyltransferase 2",
  "Intercellular adhesion molecule 1",
  "Integral membrane protein 2B",
  "Prothrombin",
  "Fibronectin",
  "Phosphatidate cytidylyltransferase 2",
  "Phosphatidate cytidylyltransferase 1",
  "Caprin-1",
  "C-type lectin domain family 10 member A",
  "cAMP-responsive element modulator",
  "Mitochondrial pyruvate carrier 2",
  "Cyclic AMP-responsive element-binding protein 3-like protein 2",
  "Kinesin-like protein KIF3B",
  "POZ-, AT hook-, and zinc finger-containing protein 1",
  "Mitochondrial pyruvate carrier 1-like protein",
  "Succinate dehydrogenase [ubiquinone] cytochrome b small subunit, mitochondrial",
  "Succinate dehydrogenase [ubiquinone] flavoprotein subunit, mitochondrial",
  "Succinate dehydrogenase [ubiquinone] iron-sulfur subunit, mitochondrial",
  "Lysosomal-associated transmembrane protein 4B",
  "IGL@ protein",
  "Succinate dehydrogenase cytochrome b560 subunit, mitochondrial",
  "Sperm-associated antigen 5",
  "Mitochondrial pyruvate carrier 1",
  "Aquaporin-5"
)

parasite_list6 <- c(
  "Integumentary mucin C.1"
)

df$broad_origin[df$group %in% human_list6]    <- "human"
df$broad_origin[df$group %in% parasite_list6] <- "parasite"

# --- Block 7 --------------------------------------------------------

human_list7 <- c(
  "Aquaporin-5",
  "Malate dehydrogenase, cytoplasmic",
  "Cathepsin S",
  "Cyclic AMP-dependent transcription factor ATF-1",
  "Fatty acid-binding protein, liver",
  "Lysosomal-associated transmembrane protein 4A",
  "Alpha-globin transcription factor CP2",
  "Aspartate aminotransferase, cytoplasmic",
  "Chymotrypsin-like elastase family member 3B",
  "Zinc finger protein 444",
  "Mannosyl-oligosaccharide 1,2-alpha-mannosidase IB",
  "Mannosyl-oligosaccharide 1,2-alpha-mannosidase IC",
  "Protein scribble homolog",
  "Nuclear factor of activated T-cells, cytoplasmic 2",
  "Fractalkine",
  "Alpha-(1,6)-fucosyltransferase",
  "Leucine-rich repeat-containing protein 7",
  "Glutamine--fructose-6-phosphate aminotransferase [isomerizing] 1",
  "DNA/RNA-binding protein KIN17",
  "Glycophorin-C",
  "Glycophorin-E",
  "V-type proton ATPase subunit d 2",
  "Pre-B-cell leukemia transcription factor 1",
  "D-3-phosphoglycerate dehydrogenase",
  "Transcriptional enhancer factor TEF-1",
  "Transcriptional enhancer factor TEF-3",
  "Transcriptional enhancer factor TEF-4",
  "Transcriptional enhancer factor TEF-5",
  "Golgin subfamily A member 5",
  "V-type proton ATPase subunit d 1",
  "Fatty acid-binding protein 5",
  "Acid sphingomyelinase-like phosphodiesterase 3b",
  "Glutamine--fructose-6-phosphate aminotransferase [isomerizing] 2",
  "Dynactin subunit 1",
  "DNA topoisomerase 3-beta-1",
  "Glycogen phosphorylase, brain form",
  "Dynamin-2",
  "Glycogen phosphorylase, liver form",
  "DNA mismatch repair protein Msh6",
  "DNA mismatch repair protein Msh2",
  "Glycogen phosphorylase, muscle form",
  "DNA mismatch repair protein Mlh1",
  "Bone marrow stromal antigen 2",
  "Mismatch repair endonuclease PMS2",
  "Protein FEV",
  "Fatty acid-binding protein, intestinal",
  "ER degradation-enhancing alpha-mannosidase-like protein 3",
  "Lamina-associated polypeptide 2, isoform alpha",
  "Lamina-associated polypeptide 2, isoforms beta/gamma",
  "Malate dehydrogenase, mitochondrial",
  "Importin subunit beta-1",
  "Fatty acid-binding protein, heart",
  "Galactosylgalactosylxylosylprotein 3-beta-glucuronosyltransferase 3",
  "Phostensin",
  "Putative aspartate aminotransferase, cytoplasmic 2",
  "Gastrotropin",
  "Amine oxidase [copper-containing] 3",
  "Amine oxidase [copper-containing] 2",
  "Uroplakin-2",
  "Fatty acid-binding protein, brain",
  "Beta-ureidopropionase",
  "Bcl-2-like protein 1",
  "Beta-secretase 1",
  "T-cell surface protein tactile",
  "Sclerostin",
  "Complement C5",
  "Mucin-13",
  "Mucin-4",
  "Plexin-B2",
  "Stromelysin-1",
  "Squalene synthase",
  "Homeobox protein TGIF2LX",
  "E3 ubiquitin ligase TRAF3IP2",
  "Cytokine SCM-1 beta",
  "Humanin-like 1",
  "Indoleamine 2,3-dioxygenase 2",
  "Indoleamine 2,3-dioxygenase 1",
  "Krueppel-like factor 10",
  "Leucine--tRNA ligase, cytoplasmic",
  "TSC22 domain family protein 1",
  "Vascular cell adhesion molecule",
  "Trefoil factor 1",
  "Neurolysin, mitochondrial",
  "Krueppel-like factor 11",
  "Lupus La protein",
  "Lymphotactin",
  "PDZ domain-containing protein 11",
  "Transcriptional regulator ERG",
  "1-phosphatidylinositol 4,5-bisphosphate phosphodiesterase gamma-2",
  "Cadherin-3",
  "Melanocyte protein PMEL",
  "Amyloid beta precursor protein binding family B member 1"
)

bacterial_list7 <- c(
  "Lipoprotein GNA1870",
  "pH-gated potassium channel KcsA"
)

toxin_list7 <- c(
  "Pesticidal crystal protein Cry2Aa"
)

df$broad_origin[df$group %in% human_list7]     <- "human"
df$broad_origin[df$group %in% bacterial_list7] <- "bacterial"
df$broad_origin[df$group %in% toxin_list7]     <- "toxin"

# --- Block 8 (final remaining 'other') ------------------------------

human_list8 <- c(
  "Amyloid beta precursor protein binding family B member 1",
  "Sulfotransferase 1E1",
  "Transmembrane 4 L6 family member 5",
  "Caspase-9",
  "Pituitary adenylate cyclase-activating polypeptide",
  "Semaphorin-4D",
  "Ectodysplasin-A",
  "Decorin",
  "Fibronectin type III domain-containing protein 1",
  "Methionine-R-sulfoxide reductase B1",
  "Methionine-R-sulfoxide reductase B2, mitochondrial",
  "Methionine-R-sulfoxide reductase B3",
  "Homeobox protein TGIF2",
  "Glycine--tRNA ligase",
  "Heat shock 70 kDa protein 1A",
  "Carboxypeptidase A4",
  "eda, eda1, eda2, ectodysplasin-a, ectodermal dysplasia protein",
  "Mitochondrial import inner membrane translocase subunit Tim17-B",
  "Mitochondrial peptide methionine sulfoxide reductase",
  "NLR family CARD domain-containing protein 4",
  "5,6-dihydroxyindole-2-carboxylic acid oxidase",
  "RNA-binding protein EWS",
  "tRNA (cytosine(34)-C(5))-methyltransferase, mitochondrial",
  "Tensin-4",
  "Prohibitin-2",
  "Hematopoietic cell signal transducer",
  "Cancer/testis antigen 1",
  "Citrate synthase, mitochondrial",
  "40-kDa huntingtin-associated protein",
  "Protein S100-A12",
  "Zinc finger C3HC-type protein 1",
  "Protein cereblon",
  "Neurotensin/neuromedin N",
  "Triosephosphate isomerase",
  "E3 ubiquitin-protein ligase RNF213",
  "Apoptosis-associated speck-like protein containing a CARD",
  "Proto-oncogene vav",
  "DNA-binding death effector domain-containing protein 2",
  "Small ribosomal subunit protein eS1",
  "Intercellular adhesion molecule 5",
  "Rho guanine nucleotide exchange factor 5",
  "X-box-binding protein 1",
  "Intercellular adhesion molecule 4",
  "Intercellular adhesion molecule 3",
  "Single-strand selective monofunctional uracil DNA glycosylase",
  "Intercellular adhesion molecule 2",
  "NACHT, LRR and PYD domains-containing protein 3",
  "DNA oxidative demethylase ALKBH2",
  "RING1 and YY1-binding protein",
  "Meiotic recombination protein DMC1/LIM15 homolog",
  "Interferon-inducible protein AIM2",
  "Cyclic AMP-responsive element-binding protein 1",
  "ATP-dependent DNA helicase Q4",
  "Lipid transferase CIDEA",
  "Lipid transferase CIDEC",
  "Plexin-B1",
  "Gelsolin",
  "Metalloproteinase inhibitor 1"
)

bacterial_list8 <- c(
  "Erythronolide synthase EryA1"
)

viral_list8 <- c(
  "Entry-fusion complex associated protein OPG095"
)

df$broad_origin[df$group %in% human_list8]     <- "human"
df$broad_origin[df$group %in% bacterial_list8] <- "bacterial"
df$broad_origin[df$group %in% viral_list8]     <- "viral"

# ---------------------- SUMMARY AFTER MANUAL OVERRIDES --------------

cat("\nUnique groups per class AFTER manual overrides (before final small block):\n")
print(table(df$broad_origin))

other_n <- sum(df$broad_origin == "other")
cat("\nRemaining 'other' groups in df (before final small block):", other_n, "\n")
if (other_n > 0) {
  cat("First remaining 'other' groups:\n")
  print(df %>% filter(broad_origin == "other") %>% head(50))
}

# --- Block 9: final manual overrides for remaining named groups ---

human_list9 <- c(
  "Glycosylphosphatidylinositol-anchored high density lipoprotein-binding protein"
)

bacterial_list9 <- c(
  "Outer surface protein A"
)

parasite_list9 <- c(
  "Polycalin",
  "Caspase Dronc",
  "Putative inactive caspase B"
)

df$broad_origin[df$group %in% human_list9]     <- "human"
df$broad_origin[df$group %in% bacterial_list9] <- "bacterial"
df$broad_origin[df$group %in% parasite_list9]  <- "parasite"


# ============================================================
# MS-RELEVANT ANTIGENS TAG (relevance_MS)
# ============================================================

# ============================================================
# MS-RELEVANT ANTIGENS TAG (relevance_MS)
# ============================================================

# Name-based lists for EBV / CMV / Akkermansia / human brain

ebv_ms_names <- c(
  "Epstein-Barr nuclear antigen 1",
  "Epstein-Barr nuclear antigen 2",
  "Epstein-Barr nuclear antigen 6",
  "Epstein-Barr nuclear antigen leader protein",
  "Envelope glycoprotein GP350",
  "G-protein coupled receptor BULF1"
)

cmv_ms_names <- c(
  "Immediate early protein IE1",
  "Serine/threonine protein kinase UL97"
)

akk_ms_names <- c(
  "Amuc_1100"
)

human_brain_names <- c(
  "Neural cell adhesion molecule L1",
  "Neural cell adhesion molecule 1",
  "Microtubule-associated protein tau",
  "Myelin-oligodendrocyte glycoprotein",
  "Sphingomyelinase",
  "Acid sphingomyelinase-like phosphodiesterase 3b"
)

# EBV / CMV / Akkermansia regex (pattern-based)
rx_ebv_ms <- regex(
  "epstein[- ]?barr|\\bebv\\b|gp350|balf1|bulf1",
  ignore_case = TRUE
)

rx_cmv_ms <- regex(
  "cytomegalovirus|\\bhcmv\\b|\\bcmv\\b|immediate early protein ie1|ul97",
  ignore_case = TRUE
)

rx_akk_ms <- regex(
  "akkermansia|amuc_?1100|\\bamuc\\b",
  ignore_case = TRUE
)

# Human brain-related antigens (heuristic list + explicit names above)
brain_terms <- c(
  "myelin basic protein",
  "myelin oligodendrocyte glycoprotein",
  "\\bmog\\b",
  "proteolipid protein 1",
  "\\bplp1\\b",
  "neurofascin",
  "contactin",
  "neurofilament",
  "glial fibrillary acidic protein",
  "\\bgfap\\b",
  "synuclein",
  "huntingtin",
  "neural cell adhesion molecule",
  "nmda receptor",
  "ampa receptor",
  "gaba receptor",
  # extra terms explicitly requested
  "neural cell adhesion molecule l1",
  "microtubule-associated protein tau",
  "sphingomyelinase",
  "acid sphingomyelinase-like phosphodiesterase 3b"
)

rx_brain_ms <- regex(paste(brain_terms, collapse = "|"), ignore_case = TRUE)

df <- df %>%
  mutate(
    relevance_MS = case_when(
      # EBV exact names or EBV-like patterns
      group %in% ebv_ms_names | str_detect(group_lc, rx_ebv_ms) ~ "EBV",
      # CMV exact names or CMV-like patterns
      group %in% cmv_ms_names | str_detect(group_lc, rx_cmv_ms) ~ "CMV",
      # Akkermansia (Amuc_1100 and similar)
      group %in% akk_ms_names | str_detect(group_lc, rx_akk_ms) ~ "Akkermansia",
      # Human brain-related antigens
      broad_origin == "human" &
        (group %in% human_brain_names | str_detect(group_lc, rx_brain_ms)) ~ "human_brain",
      TRUE ~ "none"
    )
  )

cat("\nRelevance_MS counts (group level):\n")
print(table(df$relevance_MS))

# ---------------------- JOIN BACK TO sig_ms -------------------------

sig_ms <- sig_ms %>%
  mutate(group_chr = as.character(group)) %>%
  left_join(
    select(df, group, broad_origin, relevance_MS),
    by = c("group_chr" = "group")
  ) %>%
  select(-group_chr) %>%
  mutate(label = broad_origin)

cat("\nCounts on full sig_ms rows (after overrides):\n")
print(
  sig_ms %>%
    count(label, name = "n_sequences") %>%
    arrange(desc(n_sequences))
)

cat("\nRelevance_MS counts on full sig_ms:\n")
print(
  sig_ms %>%
    count(relevance_MS, name = "n_sequences") %>%
    arrange(desc(n_sequences))
)


# ---------------------- WRITE OUTPUT --------------------------------

write_tsv(sig_ms, out_path)
cat("\nWritten labeled results to:", out_path, "\n")

