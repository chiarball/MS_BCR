library(dplyr, lib.loc = "/doctorai/niccoloc/libR2")
library(stringr, lib.loc = "/doctorai/niccoloc/libR2")
library(purrr, lib.loc = "/doctorai/niccoloc/libR2")
library(tidyr , lib.loc = "/doctorai/niccoloc/libR2")
library(stringdist , lib.loc = "/doctorai/niccoloc/libR2")
library(viridis , lib.loc = "/doctorai/niccoloc/libR2")
library(isoband, lib.loc = "/doctorai/niccoloc/libR2")
library(data.table, lib.loc = "/doctorai/niccoloc/libR2")

# pgen <- fread("/doctorai/mariaabb/PER_CHIARA_df_sharing_malati_withpgen.tsv")
pgen <- fread("/doctorai/niccoloc/MS_db/MS_BCR/df_pgen.tsv")%>% 
  select(-tretment,-cell_type)

pgen %>% 
  mutate(clone_class=case_when(
    # top_cluster==TRUE ~ "TOP_Cluster",
    presence_group=="MS_only" & MS_only_cluster==TRUE ~ "MS_cluster",
    presence_group=="HC_only"  & HC_only_cluster==T & is_public==T~ "HC_cluster",
    presence_group=="MS_only" & MS_only_cluster==FALSE ~ "MS_public",
    presence_group=="MS_and_HC" ~ "public",
    presence_group=="HC_only"  & is_public==T~ "HC_public",
    .default =   "private"
    
  )) -> pgen

table(pgen$clone_class)

# Step 1: filtro solo le junction_aa NON presenti nei sani
MS_only <- pgen %>%
  filter(clone_class == "MS_cluster",
         Pgen != 0)

MS_HC_public <- pgen %>%
  filter(clone_class == "public",
         Pgen != 0)



pgen_thr=1e-16

# Step 2: conto quante hanno Pgen < 1e-10
n_low_pgen <- MS_only %>%
  filter(Pgen < pgen_thr) %>%
  nrow()

# Step 3: totale delle junction_aa non presenti nei sani
total_false <- nrow(MS_only)

# Step 4: proporzione
prop_low_pgen <- n_low_pgen / total_false

cat("Number with Pgen < ",pgen_thr, n_low_pgen, "\n")
cat("Total junction_aa with present_in_sani == FALSE:", total_false, "\n")
cat("Proportion:", prop_low_pgen, "\n")


#present in sani 




### plotting
install.packages('kSamples', lib="/doctorai/niccoloc/libR2")
library(kSamples,lib.loc="/doctorai/niccoloc/libR2")
library(ggpubr, lib.loc="/doctorai/chiarba/lib")


p1=pgen %>% 
  filter(clone_class %in% c('MS_cluster','MS_and_HC', 'public')) %>%
  mutate(Pgen=ifelse(Pgen==0, 1e-50, Pgen))

 table(p1$clone_class, p1$is_zero) 

ad.test(
  x = p1$Pgen[p1$clone_class == "HC_cluster"],
  y = p1$Pgen[p1$clone_class == "MS_cluster"]
)

by(
  log10(p1$Pgen),
  p1$clone_class,
  summary
)



# Stastistical test: Wilcoxon  
p1=pgen %>% 
  filter(clone_class %in% c('MS_cluster','public')) %>%
  mutate(Pgen=ifelse(Pgen==0, 1e-50, Pgen))
table(p1$clone_class)


p_test <- wilcox.test(Pgen ~ clone_class, data = p1)

p_value <- p_test$p.value
p_value


p2=ggplot(p1, aes(x = Pgen, color = clone_class)) +
  geom_density(   ) +
  # geom_histogram( aes(y=..density.., fill=clone_class), position="identity", alpha=0.4, bins=50 ) +
  scale_x_log10( 
    limits = c(1e-35, NA),
    breaks= c(1e-30, 1e-25,1e-20,1e-15, 1e-10, 1e-5 )
                 ) +
  # scale_fill_manual(values = c("steelblue", "orange")) +
  scale_color_manual(values = c("orange", "steelblue"), labels = c("MS cluster only", "MS and HC")) +
  labs(
    x = "Pgen (log10 scale)",
    y = "Density",
    color = "Present in:"
  ) +
  
  theme_bw(base_size = 14)
p2

### stat wilcox test

wilcox_result <- ks.test(MS_only$Pgen,
                         MS_HC_public$Pgen )

wilcox_result <- wilcox.test(MS_only$Pgen,
                             MS_HC_public$Pgen,
                             alternative =  'less')


p_value <- wilcox_result$p.value
W_stat  <- wilcox_result$statistic

median_MS_only  <- median(MS_only$Pgen)
median_MS_HC    <- median(MS_HC_public$Pgen)

iqr_MS_only     <- IQR(MS_only$Pgen)
iqr_MS_HC       <- IQR(MS_HC_public$Pgen)

print(wilcox_result)
print(paste("p-value:", p_value))
print(paste("W statistic:", W_stat))
print(paste("Median MS-only:", median_MS_only))
print(paste("Median MS+HC:", median_MS_HC))
print(paste("IQR MS-only:", iqr_MS_only))
print(paste("IQR MS+HC:", iqr_MS_HC))

grDevices::svg(
  filename = "/doctorai/niccoloc/pgen_MS_db.svg",
  width    = 9,
  height   = 7
)
print(p2)
dev.off()


