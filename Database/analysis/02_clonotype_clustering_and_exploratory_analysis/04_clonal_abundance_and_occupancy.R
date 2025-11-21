## 04_clonal_abundance_and_occupancy.R
## Rare clonal proportion and cumulative Top-N clonotypes + quantitative summaries

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(tibble)
library(purrr)

## ------------------------------------------------------------------
## 1) Read data and build abundance bins (bin_freq)
## ------------------------------------------------------------------

# Read per-patient clustered table
df <- fread("/doctorai/chiarba/analysis/MS_HC_naive_per_patient_h0h1.tsv")

# Harmonize tissue label for HC blood
df[tissue == "Peripheral blood", tissue := "blood from HC"]

unique(df$tissue)

# Count clones per subject × tissue × cluster_h1
clone_counts <- df %>%
  group_by(subject_number, tissue, cluster_h1) %>%
  summarise(count = n(), .groups = "drop")

# Assign each clone to a “count-bin”
clone_counts <- clone_counts %>%
  mutate(count_bin = cut(
    count,
    breaks = c(1, 2, 4, 11, 31, 101, Inf),
    labels = c("1", "2–3", "4–10", "11–30", "31–100", "101+"),
    right  = FALSE
  ))

# Compute relative proportion per subject × tissue × bin
bin_freq <- clone_counts %>%
  group_by(subject_number, tissue, count_bin) %>%
  summarise(bin_sum  = sum(count), .groups = "drop") %>%
  group_by(subject_number, tissue) %>%
  mutate(rel_prop = bin_sum / sum(bin_sum)) %>%
  ungroup()

# ensure tissue is a character vector (avoid accidental factor issues)
bin_freq$tissue <- as.character(bin_freq$tissue)

# show actual unique tissue names so you can spot typos
unique_tissues <- sort(unique(bin_freq$tissue))
print(unique_tissues)

# fixed first facets you want on top, in that order (adjust text if needed)
fixed_first <- c("blood form HC", "blood", "CLN", "CSF", "Pia mater",
                 "Choroid plexus", "Brain lesion")

# approximate-match each requested name against actual tissue names
matched <- sapply(
  fixed_first,
  function(x) {
    # exact match preferred
    if (x %in% unique_tissues) return(x)
    # approximate match (first match)
    approx <- unique_tissues[agrep(x, unique_tissues, max.distance = 0.2)]
    if (length(approx) >= 1) return(approx[1])
    NA_character_
  },
  USE.NAMES = FALSE
)

# display mapping
mapping <- data.frame(
  requested = fixed_first,
  matched   = matched,
  stringsAsFactors = FALSE
)
print(mapping)

# if any requested names failed to match, warn and drop NAs
if (any(is.na(matched))) {
  warning("Some requested facet names were not found exactly. Check 'mapping' and adjust 'fixed_first' if needed.")
  matched <- matched[!is.na(matched)]
}

# compute totals per tissue
totals <- aggregate(rel_prop ~ tissue, data = bin_freq, FUN = sum, na.rm = TRUE)

# exclude the already-fixed ones and order the rest by descending total rel_prop
rest <- totals[!totals$tissue %in% matched, , drop = FALSE]
rest_ordered <- as.character(rest[order(-rest$rel_prop), "tissue"])

# final levels: fixed first (in the requested order), then the ordered rest
new_levels <- c(matched, rest_ordered)

# apply to the dataframe
bin_freq$tissue <- factor(bin_freq$tissue, levels = new_levels)

# store desired_levels so that the Top-N plot can reuse the same facet order
desired_levels <- levels(bin_freq$tissue)

# Plot 1: stacked bars of abundance bins per patient and tissue
ggplot(bin_freq,
       aes(x = factor(subject_number),
           y = rel_prop,
           fill = count_bin)) +
  geom_col(color = "black", size = 0.2) +
  facet_wrap(~ tissue, scales = "free_x", ncol = 4) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_brewer("Clonotype count", palette = "Spectral") +
  labs(x = "Patient", y = "Occupied repertoire (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        legend.title = element_text(face = "bold"))

## ------------------------------------------------------------------
## 2) Cumulative Top-N clonotypes (requires clone_per_clone)
## ------------------------------------------------------------------
## Assumes you already have `clone_per_clone` with columns:
##   subject_number, tissue, freq  (freq = relative frequency of each clone)

topN <- c(1, 5, 10, 20, 50)

# helper: lightens colors by mixing them with white
lighten <- function(cols, frac = 0.4) {
  # cols: vector of hex colors; frac: 0..1 fraction of white to mix in
  rgb_mat <- grDevices::col2rgb(cols)
  white <- 255
  mixed <- round((1 - frac) * rgb_mat + frac * white)
  grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ], maxColorValue = 255)
}

# base palette (keeps blues) and then lighten
pal_base <- brewer.pal(length(topN) + 1, "Spectral")
pal <- lighten(pal_base, frac = 0.45)

# compute cumulative per clone
cum <- clone_per_clone %>%
  group_by(subject_number, tissue) %>%
  arrange(desc(freq), .by_group = TRUE) %>%
  mutate(
    cumfreq = cumsum(freq),
    rank    = row_number()
  ) %>%
  ungroup()

# for each sample compute the incremental parts: Top1, (2-5), (6-10), (11-20), (21-50), >50
parts <- cum %>%
  group_by(subject_number, tissue) %>%
  group_modify(~ {
    df_sub <- .x
    maxcum <- if (nrow(df_sub) > 0) max(df_sub$cumfreq, na.rm = TRUE) else 0
    # cumulative at requested thresholds: if missing use maxcum
    cum_at <- sapply(topN, function(n) {
      v <- df_sub$cumfreq[df_sub$rank == n]
      if (length(v) == 0) maxcum else v
    })
    # incremental parts
    inc <- c(cum_at[1], diff(cum_at), pmax(0, 1 - cum_at[length(cum_at)]))
    tibble(N = c(paste0("Top ", topN), ">50"), value = inc)
  }, .keep = TRUE) %>%
  ungroup()

# ensure factor levels stable and tissue ordering preserved
parts <- parts %>%
  mutate(
    N      = factor(N, levels = c(paste0("Top ", topN), ">50")),
    tissue = factor(tissue, levels = desired_levels)
  )

# Plot 2: stacked bars summing to 100% (Top-N cumulative)
ggplot(parts, aes(x = factor(subject_number), y = value, fill = N)) +
  geom_col(color = "black", size = 0.2) +
  facet_wrap(~ tissue, scales = "free_x", ncol = 4) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual("Cumulative", values = pal) +
  labs(
    x = "Sample",
    y = "Cumulative proportion",
    title = "Cumulative repertoire captured by Top-N clonotypes",
    subtitle = "Segments: Top1, Top2–5, Top6–10, Top11–20, Top21–50, >50"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        legend.title = element_text(face = "bold"))

## ------------------------------------------------------------------
## 3) Quantitative summaries MS vs HC, CSF vs blood
## ------------------------------------------------------------------

##### 3a) Abundance-binned summary (clone count bins) #####

# bin_freq_q: add cohort and tissue_group (CSF vs blood)
bin_freq_q <- bin_freq %>%
  mutate(
    cohort = if_else(grepl("^HC", subject_number), "HC", "MS"),
    tissue_group = case_when(
      grepl("csf", tissue,   ignore.case = TRUE) ~ "CSF",
      grepl("blood", tissue, ignore.case = TRUE) ~ "blood",
      TRUE                                       ~ NA_character_
    )
  ) %>%
  filter(!is.na(tissue_group))

# parts_q: same idea (cohort + tissue_group)
parts_q <- parts %>%
  mutate(
    cohort = if_else(grepl("^HC", subject_number), "HC", "MS"),
    tissue_group = case_when(
      grepl("csf", tissue,   ignore.case = TRUE) ~ "CSF",
      grepl("blood", tissue, ignore.case = TRUE) ~ "blood",
      TRUE                                       ~ NA_character_
    )
  ) %>%
  filter(!is.na(tissue_group))

# For each cohort (MS/HC), tissue_group (CSF/blood), and count_bin:
# - mean and median proportion across patients
# - number of patients in that group
bin_summary <- bin_freq_q %>%
  group_by(cohort, tissue_group, count_bin) %>%
  summarise(
    mean_prop   = mean(rel_prop),
    median_prop = median(rel_prop),
    sd_prop     = sd(rel_prop),
    n_patients  = n_distinct(subject_number),
    .groups     = "drop"
  ) %>%
  mutate(
    mean_pct   = round(mean_prop   * 100, 1),
    median_pct = round(median_prop * 100, 1)
  ) %>%
  arrange(cohort, tissue_group, count_bin)

bin_summary

##### 3b) Rank-based cumulative Top-N summary #####

# For each cohort (MS/HC), tissue_group (CSF/blood), and N ("Top 1", "Top 2-5", ..., ">50"):
# - mean and median fraction across patients
# - number of patients in that group

topN_summary <- parts_q %>%
  group_by(cohort, tissue_group, N) %>%
  summarise(
    mean_prop   = mean(value),
    median_prop = median(value),
    sd_prop     = sd(value),
    n_patients  = n_distinct(subject_number),
    .groups     = "drop"
  ) %>%
  mutate(
    mean_pct   = round(mean_prop   * 100, 1),
    median_pct = round(median_prop * 100, 1)
  ) %>%
  arrange(cohort, tissue_group, N)

topN_summary
