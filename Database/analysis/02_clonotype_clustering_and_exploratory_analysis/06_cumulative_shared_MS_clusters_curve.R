cluster_db <- fread("/doctorai/chiarba/analysis/MS_HC_naive_h0h1.tsv")

library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Mapping DB name 

study_map <- tibble::tribble(
  ~study_db,           ~study,
  "Palanichemy",       "Palanichamy et al.",
  "Stern",             "Stern et al.",
  "Greenfield",        "Greenfield AL",
  "Ramesh",            "Ramesh et al.",
  "lomakin",           "Lomakin et al.",
  "Agrafiotis",        "Agrafiotis et al.",
  "Laurent",           "Laurent S.A. et al.",
  "Perez_Salvidar",    "Pérez-Saldívar et al.",
  "Ruschill",          "Ruschil et al."
)


# 2) Year table 

study_year <- tibble::tribble(
  ~study,                 ~anno,
  "Stern et al.",          2014,
  "Palanichamy et al.",    2014,
  "Greenfield AL",         2019,
  "Ramesh et al.",         2020,
  "Lomakin et al.",        2022,
  "Agrafiotis et al.",     2023,
  "Laurent S.A. et al.",   2023,
  "Ruschil et al.",        2023,
  "Pérez-Saldívar et al.", 2024
)

# Order
study_order <- study_year %>% arrange(anno, study) %>% pull(study)


# 3) Base DB (include tissue) + filter cluster touching ^H

db0 <- cluster_db %>%
  select(cluster_h1, subject_number, study, tissue) %>%
  mutate(study = as.character(study))

clusters_with_H <- db0 %>%
  filter(grepl("^H", subject_number)) %>%
  distinct(cluster_h1) %>%
  pull(cluster_h1)

db <- db0 %>%
  filter(!cluster_h1 %in% clusters_with_H)

# 4) Normalize study names + keep only mapped/present

db <- db %>%
  rename(study_db = study) %>%
  left_join(study_map, by = "study_db") %>%
  filter(!is.na(study)) %>%
  mutate(study = factor(study, levels = study_order))

studies <- levels(db$study)


# 5) Build axis labels: Study + Year + Tissues

tissue_by_study <- db %>%
  distinct(study, tissue) %>%
  group_by(study) %>%
  summarise(
    tissues = {
      tt <- sort(unique(tissue))
      if (length(tt) <= 2) {
        paste(tt, collapse = " / ")
      } else {
        n <- length(tt)
        top <- tt[1:ceiling(n / 2)]
        bottom <- tt[(ceiling(n / 2) + 1):n]
        paste0(
          paste(top, collapse = " / "),
          "\n",
          paste(bottom, collapse = " / ")
        )
      }
    },
    .groups = "drop"
  )


axis_labels <- study_year %>%
  mutate(study = factor(study, levels = study_order)) %>%
  left_join(tissue_by_study, by = "study") %>%
  mutate(
    tissues = replace_na(tissues, "NA"),
    xlab = paste0(as.character(study), "\n", anno, "\n", tissues)
  ) %>%
  arrange(study)

label_vec <- setNames(axis_labels$xlab, as.character(axis_labels$study))


# 6) WITHIN (per-studio): # cluster_h1 shared by >=2 subjects within the same study

within_by_study <- db %>%
  group_by(study, cluster_h1) %>%
  summarise(n_subj = n_distinct(subject_number), .groups = "drop") %>%
  filter(n_subj >= 2) %>%
  count(study, name = "within_clusters")


# 7) ACROSS (cumulative): # cluster_h1 present in >=2 studies in the cumulative set

across_list <- lapply(seq_along(studies), function(i) {
  db_i <- db %>% filter(study %in% studies[1:i])

  n_clusters_i <- db_i %>%
    group_by(cluster_h1) %>%
    summarise(n_study = n_distinct(study), .groups = "drop") %>%
    filter(n_study >= 2) %>%
    nrow()

  tibble(
    study = factor(studies[i], levels = studies),
    across_clusters_cum = n_clusters_i
  )
})

across_cum <- bind_rows(across_list)


# 7b) TOTAL (cumulative): union of within + across, no double-count
#     = clusters with >=2 subjects in the cumulative set (anywhere)

total_list <- lapply(seq_along(studies), function(i) {
  db_i <- db %>% filter(study %in% studies[1:i])

  n_total_i <- db_i %>%
    group_by(cluster_h1) %>%
    summarise(n_subj = n_distinct(subject_number), .groups = "drop") %>%
    filter(n_subj >= 2) %>%
    nrow()

  tibble(
    study = factor(studies[i], levels = studies),
    total_shared_clusters_cum = n_total_i
  )
})

total_cum <- bind_rows(total_list)


# Number of sequences per study (post ^H filter)


study_sizes <- db %>%
  count(study, name = "n_sequences_noH") %>%
  mutate(
    study_size = paste0(
      as.character(study),
      " (",
      format(n_sequences_noH, big.mark = ",", trim = TRUE),  
      ")"
    )
  )


# 8) Plot table + plot

plot_df <- tibble(study = factor(studies, levels = studies)) %>%
  left_join(within_by_study, by = "study") %>%
  left_join(across_cum, by = "study") %>%
  left_join(total_cum, by = "study") %>%
  left_join(
    study_sizes %>% select(study, study_size),
    by = "study"
  ) %>%
  mutate(
    within_clusters = replace_na(within_clusters, 0L),
    across_clusters_cum = replace_na(across_clusters_cum, 0L),
    total_shared_clusters_cum = replace_na(total_shared_clusters_cum, 0L)
  ) %>%
  pivot_longer(
    cols = c(
      within_clusters,
      across_clusters_cum,
      total_shared_clusters_cum
    ),
    names_to = "linea",
    values_to = "n_clusters"
  )

size_values <- c(
  "Agrafiotis et al. (880)"            = 1.5,
  "Lomakin et al. (35,538)"            = 2.0,
  "Laurent S.A. et al. (54,344)"       = 2.3,
  "Ramesh et al. (87,907)"             = 2.7,
  "Palanichamy et al. (120,207)"       = 3.0,
  "Pérez-Saldívar et al. (125,313)"    = 3.2,
  "Ruschil et al. (155,954)"           = 3.5,
  "Greenfield AL (409,405)"            = 4.5,
  "Stern et al. (539,607)"             = 5
)

ggplot(
  plot_df,
  aes(
    x = study,
    y = n_clusters,
    group = linea,
    color = linea,
    size = study_size
  )
) +
  geom_line(size = 1) +
  geom_point(alpha = 0.75) +
  scale_x_discrete(labels = label_vec) +
  scale_y_continuous(
  limits = c(0, NA),
  expand = expansion(mult = c(0.05, 0.05))  
)+
  scale_color_manual(
    values = c(
      "within_clusters"           = "#0072B2",
      "across_clusters_cum"       = "#E69F00",
      "total_shared_clusters_cum" = "#009E73"
    ),
    labels = c(
      "within_clusters"           = "within study",
      "across_clusters_cum"       = "across studies",
      "total_shared_clusters_cum" = "total shared"
    ),
    guide = guide_legend(
      ncol = 1,        
      byrow = TRUE,
      override.aes = list(size = 3)
    )
  ) +
scale_size_manual(
  values = size_values,
  name = "Sequences per study",
  guide = guide_legend(
    ncol = 1,
    byrow = TRUE,
    override.aes = list(alpha = 1)
  )
) +
  labs(
    x = "Study/Year/Tissues",
    y = "Number of shared cluster only in MS patients",
    color = NULL
  ) +
  theme_bw() +
  theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "right",     
  legend.box = "vertical",       
  legend.spacing.y = unit(0.3, "cm"),
  legend.key.width = unit(0.9, "cm")
)

