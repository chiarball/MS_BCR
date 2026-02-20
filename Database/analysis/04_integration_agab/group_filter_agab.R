#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr, lib.loc = "/doctorai/niccoloc/libR2")
  library(stringr, lib.loc = "/doctorai/niccoloc/libR2")
  library(purrr, lib.loc = "/doctorai/niccoloc/libR2")
  library(tidyr,  lib.loc = "/doctorai/niccoloc/libR2")
  library(stringdist, lib.loc = "/doctorai/niccoloc/libR2")
  library(tidytext, lib.loc = "/doctorai/niccoloc/libR2")
  library(data.table, lib.loc = "/doctorai/niccoloc/libR2")
  library(withr, lib.loc = "/doctorai/niccoloc/libR2")
  library( arrow , lib = "/doctorai/niccoloc/libR2"   )
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: score_x2_group.R <input_path> <output_path> <group_by_vgene>",
    call. = FALSE
  )
}
print("Arguments:")
print(args)

input_path  <- args[1]
output_path <- args[2]
vgene_group <- ifelse(length(args) >= 3, as.logical(args[3]), FALSE)
# cluster_h   <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.20)

print(paste("Input path:", input_path))
print(paste("Output path:", output_path))

# -----------------------------
# 1) Helpers: parsing + cleaning
# -----------------------------

# Split "a; b; c" into character vector, trimmed, unique, drop empties
split_antigens <- function(which_ag, sep = ";") {
  if (is.na(which_ag) || which_ag == "") return(character(0))
  which_ag %>%
    str_split(paste0("\\s*", sep, "\\s*")) %>%
    .[[1]] %>%
    str_trim() %>%
    discard(~ .x == "" || is.na(.x)) %>%
    unique()
}

# A simple but effective blacklist for "not real antigen" / scaffold noise
# Tune as you like
is_blacklisted <- function(x) {
  pat <- paste0(
    "(",
    paste(c(
      "immunoglobulin",
      "igg\\b", "ige\\b", "iga\\b", "igm\\b", "igd\\b",
      "heavy constant", "light constant",
      "\\bfc\\b", "fragment crystallizable", "antibody",
      "constant region", "constant gamma", "constant epsilon", "constant alpha",
      "kappa constant", "lambda constant",
      "t cell receptor beta constant", "t cell receptor alpha constant"
    ), collapse = "|"),
    ")"
  )
  str_detect(str_to_lower(x), pat)
}

# Normalize antigen strings to reduce typos/format differences.
# Keeps it conservative (won't over-collapse unrelated stuff).
normalize_ag <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("\\(.*?\\)", " ") %>%                 # remove parentheses content
    str_replace_all("[^a-z0-9\\s\\-_/]+", " ") %>%        # remove punctuation (keep - _ /)
    str_replace_all("\\b(human|mouse|murine|rat)\\b", " ") %>%
    str_replace_all("\\b(protein|receptor|antigen|precursor|molecule)\\b", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

# Optional alias normalization for common immune checkpoint / frequent targets.
# Expand as needed for your domain.
canonicalize_aliases <- function(x) {
  # order matters
  x %>%
    str_replace_all("\\bpd\\-?1\\b|\\bpdcd1\\b|\\bcd279\\b", "pd1") %>%
    str_replace_all("\\bpd\\-?l1\\b|\\bcd274\\b|\\bb7\\-?h1\\b", "pdl1") %>%
    str_replace_all("\\bpd\\-?l2\\b|\\bpdcd1lg2\\b|\\bcd273\\b", "pdl2") %>%
    str_replace_all("\\bctla\\-?4\\b|\\bcd152\\b", "ctla4") %>%
    str_replace_all("\\bher\\-?2\\b|\\berbb2\\b", "her2") %>%
    str_replace_all("\\bher\\-?1\\b|\\begfr\\b", "egfr") %>%
    str_replace_all("\\bvegf\\-?a\\b|\\bvegfa\\b", "vegfa") %>%
    str_replace_all("\\btnf\\-?a\\b|\\btnfa\\b", "tnfa") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

# Apply full cleaning pipeline to an antigen vector
clean_antigens <- function(ag_vec,
                           drop_blacklist = TRUE,
                           do_aliases = TRUE) {
  if (length(ag_vec) == 0) return(character(0))
  ag <- ag_vec %>% discard(~ is.na(.x) || .x == "")
  if (drop_blacklist) ag <- ag %>% discard(is_blacklisted)
  ag <- normalize_ag(ag)
  if (do_aliases) ag <- canonicalize_aliases(ag)
  ag <- ag %>% discard(~ .x == "" || is.na(.x)) %>% unique()
  ag
}



# --- helpers ---
normalize_text <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("5'\\-", "5 ") %>%     # 5'- -> "5 "
    str_replace_all("[^a-z0-9\\s]", " ") %>%
    str_squish()
}

default_stop <- c(
  "the","a","an","and","or","of","to","in","on","by",
  "cell","cells","surface","protein","glycoprotein","receptor","expressed",
  "human","mouse","bovin","rat","pig"
)

tokenize_labels <- function(labels, stop_words = default_stop) {
  tibble(label = labels, doc_id = seq_along(labels)) %>%
    mutate(text = normalize_text(label)) %>%
    unnest_tokens(token, text, token = "words") %>%
    filter(!token %in% stop_words, nchar(token) > 2)
}

bow_top_terms <- function(labels, top_k = 8, stop_words = default_stop, ngrams = 1) {
  labels <- labels[!is.na(labels) & str_trim(labels) != ""]
  if (length(labels) == 0) return(tibble(term = character(0), freq = integer(0)))
  
  df <- tibble(label = labels) %>%
    mutate(text = normalize_text(label))
  
  if (ngrams == 1) {
    out <- df %>%
      unnest_tokens(term, text, token = "words") %>%
      filter(!term %in% stop_words, nchar(term) > 2) %>%
      count(term, sort = TRUE) %>%
      slice_head(n = top_k)
  } else {
    out <- df %>%
      unnest_tokens(term, text, token = "ngrams", n = ngrams) %>%
      mutate(term = str_squish(term)) %>%
      # remove ngrams that contain stopwords (simple rule)
      filter(!str_detect(term, paste0("\\b(", paste(stop_words, collapse = "|"), ")\\b"))) %>%
      count(term, sort = TRUE) %>%
      slice_head(n = top_k)
  }
  out
}



# -----------------------------
# 2) Coherence metrics
# -----------------------------

# Lexical coherence: mean pairwise (1 - distance) using Jaro-Winkler
lexical_coherence <- function(ag_vec) {
  k <- length(ag_vec)
  if (k <= 1) return(1)
  d <- stringdistmatrix(ag_vec, ag_vec, method = "jw")
  sim <- 1 - d
  mean(sim[upper.tri(sim)], na.rm = TRUE)
}

# Cluster purity: hierarchical clustering on JW distance
# purity = size(largest cluster)/k
semantic_purity <- function(ag_vec, h = 0.20) {
  k <- length(ag_vec)
  if (k <= 1) return(1)
  d <- stringdistmatrix(ag_vec, ag_vec, method = "jw")
  hc <- hclust(as.dist(d), method = "average")
  cl <- cutree(hc, h = h)
  as.numeric(max(table(cl)) / k)
}

# Entropy-based support from frequency distribution. --- dont use it
entropy_support <- function(ag_vec) {
  # If you only have unique antigens (as in x2_group$which_ag), entropy is maximal.
  # Better version uses raw antigen counts from the original x2 rows (see section 4).
  k <- length(ag_vec)
  if (k <= 1) return(1)
  
  # With unique list: all p_i = 1/k -> H = log(k) -> normalized entropy = 1 -> support = 0
  # Still return that, but you'll likely want the frequency-aware version below.
  p <- rep(1 / k, k)
  H <- -sum(p * log(p))
  Hnorm <- H / log(k)
  1 - Hnorm
} 

# Frequency-aware entropy support: takes a vector with repeats (raw antigens, not unique) --- dont use it
entropy_support_freq <- function(raw_ag_vec) {
  raw_ag_vec <- raw_ag_vec %>% discard(~ is.na(.x) || .x == "")
  if (length(raw_ag_vec) <= 1) return(1)
  tab <- table(raw_ag_vec)
  p <- as.numeric(tab) / sum(tab)
  k <- length(p)
  if (k <= 1) return(1)
  H <- -sum(p * log(p))
  Hnorm <- H / log(k)
  1 - Hnorm
}

# Final score: weighted + evidence scaling ---- dont use it
final_score <- function(lexical, purity, support, n,
                        w1 = 0.3, w2 = 0.4, w3 = 0.3) {
  base <- w1 * lexical + w2 * purity + w3 * support
  base * log1p(n)
}

# -----------------------------
# 3) Main scorer for x2_group
# -----------------------------

score_x2_group <- function(x2_group,
                           cluster_h = 0.20,
                           drop_blacklist = TRUE,
                           do_aliases = TRUE,
                           w1 = 0.3, w2 = 0.4, w3 = 0.3) {
  
  x2_group %>%
    mutate(
      # parse which_ag
      ag_raw_unique = map(which_ag, split_antigens),
      # clean / normalize
      ag_clean = map(ag_raw_unique, ~ clean_antigens(.x,
                                                     drop_blacklist = drop_blacklist,
                                                     do_aliases = do_aliases)),
      n_antigens_clean = map_int(ag_clean, length),
      
      lexical = map_dbl(ag_clean, lexical_coherence),
      purity  = map_dbl(ag_clean, ~ semantic_purity(.x, h = cluster_h)),
      
      # With x2_group you only have unique antigens -> entropy_support will often be low.
      # Keep it for completeness; prefer the frequency-aware version in section 4.
      support = map_dbl(ag_clean, entropy_support),
      
      coherence = (w1 * lexical + w2 * purity + w3 * support),
      final_score = final_score(lexical, purity, support, n, w1, w2, w3),
      
      # rebuilt canonical antigen list (cleaned)
      which_ag_clean = map_chr(ag_clean, ~ paste(.x, collapse = "; "))
    ) %>%
    arrange(desc(final_score))
}




x2=fread(input_path)

x2=x2 %>% 
  mutate(target_name= ifelse( target_name=="" ,
                              str_c(target_pdb,target_uniprot,sep = ''),target_name)) %>%
  mutate(antigen= case_when(
    dataset=='hiv' ~ 'HIV',
    TRUE ~ antigen
  ),
  target_name= case_when(
    dataset=='hiv' ~ 'HIV',
    TRUE ~ target_name
  ))%>% distinct() 

# test2= x2 %>%
#   filter(confidence != 'medium', dataset!='hiv', target_name!="") %>% 
#   group_by(heavy_sequence,light_sequence,dataset,target_name ) %>%
#   summarise(n= n() 
#             # target_pdb,target_uniprot
#   )
# View(test2 %>% filter(n>1) %>% arrange(desc(n)))


test2= x2 %>%
  filter(confidence != 'medium', dataset!='hiv', target_name!="") %>% 
  group_by(cdr3_aa,dataset,antigen ) %>%
  summarise(n= n() )

non_patent_db_unique=x2 %>%
  filter(confidence != 'medium',
         dataset!='hiv',
         target_name!="") %>% 
  distinct(cdr3_aa,dataset,antigen , .keep_all = TRUE)
  # group_by(cdr3_aa,dataset,antigen ) %>%
  # summarise(n= n() )   %>% View()
if(vgene_group){
  x2_group= x2  %>% 
    filter(confidence == 'medium')  %>%
    mutate(mut_stat= ifelse( ighv_pident>=98,'unmutated','mutated')  ) %>%
    group_by(cdr3_aa,ighv_v_gene,mut_stat )  %>% 
    summarise(n= n( )  ,
              n_antigens= n_distinct(antigen),
              which_ag= paste0( unique(antigen), collapse = '; '  ),
              target_ag=paste0( unique(target_name), collapse = '; '  )
    )
  
  non_patent_db_unique=x2 %>%
    filter(confidence != 'medium', 
           # dataset!='hiv',
           target_name!="") %>% 
    mutate(mut_stat= ifelse( ighv_pident>=98,'unmutated','mutated')  ) %>%
    distinct(cdr3_aa,ighv_v_gene,mut_stat,dataset,antigen , .keep_all = TRUE)
  
  
  
} else{
x2_group= x2  %>% 
  #select only medium confidence, which is the patents database
  filter(confidence == 'medium')  %>% 
  group_by(cdr3_aa )  %>% 
  summarise(n= n( )  ,
            n_antigens= n_distinct(antigen),
            which_ag= paste0( unique(antigen), collapse = '; '  ),
            target_ag=paste0( unique(target_name), collapse = '; '  )
        
  )

non_patent_db_unique=x2 %>%
  filter(confidence != 'medium', dataset!='hiv', target_name!="") %>% 
  distinct(cdr3_aa,dataset,antigen , .keep_all = TRUE)


}

print("read evertyhing")
#early stopping for debug
# quit(save = "no", status = 0)


# Score and cluster by lexical JW similarity
scored <- score_x2_group(x2_group, cluster_h = 0.20, drop_blacklist = F)

#print for debugs
print("scored head:")
print(head(scored))

scored_ag_2= scored  
  # mutate(
  #   perc_ag_clean = n_antigens_clean / n_antigens * 100
  # ) 
# # %>% 
#   filter(
#     # perc_ag_clean >=80 ,
#           
#           n_antigens >1   )  


unique_ags = scored_ag_2 %>%     pull (target_ag) %>% str_split('; ') %>% unlist() %>% unique()
unique_ags_species = unique_ags%>% str_remove('^.*_') %>% unique()    

#get the unique genes without species suffix, from the provided target_ag column from the patents db authors
scored2 <- scored_ag_2 %>%
  mutate(
    target_ag_clean = map(
      str_split(target_ag, ";\\s*"),
      #remove the species suffixes from unique_ags_species vector
      ~ str_remove(.x , paste0("_(", paste(unique_ags_species, collapse = "|"), ")$") )),
    perc_identical_target = map_dbl(
      target_ag_clean,
      ~ {
        if (length(.x) == 0) return(NA_real_)
        tab <- table(.x)
        max(tab) / sum(tab) * 100
      }  )
  )
#print for debugs
print("scored2 head:")
print(head(scored2))

patents_filtered=scored2 %>% 
  mutate(patent_quality=case_when(
    n_antigens ==1 ~ 'high', # we keep all with 1 antigen
    n_antigens <=2 ~ 'good', # we keep all with 2 or less antigens
    purity >=0.7 ~ 'good', #if more than 2 antigens, we keep the groups with high lexical clustering purity
    purity <0.7 & purity >=0.5 & n_antigens <=10 ~ 'medium', # if purity is medium, only the ones with less than 10 annotated ags
    .default = 'low'
  ),
  keep=ifelse(patent_quality  == 'low',FALSE,TRUE)  ) %>%
  filter(keep) 


#use a bag of words approach to extract top terms from the antigen list
df_out <- patents_filtered %>%
  mutate(
    bow_terms = map(which_ag, bow_top_terms, top_k = 8, ngrams = 1)
  )  

#print for debugs
print("df_out head:")
print(head(df_out))

# Add "generic terms"
GENERIC <- c( "protein","molecule")


make_label_maxanchored <- function(bow_tbl,
                                   p = 0.7,      # keep terms with n >= p * max_n
                                   min_n = 2,    # also require an absolute minimum
                                   generic = GENERIC,
                                   collapse = "_") {
  if (is.null(bow_tbl) || nrow(bow_tbl) == 0) return(NA_character_)
  
  x <- bow_tbl %>%
    mutate(term = str_to_lower(term)) %>%
    filter(!term %in% generic)
  
  if (nrow(x) == 0) return(NA_character_)
  
  max_n <- max(x$n, na.rm = TRUE)
  thr <- max(min_n, ceiling(p * max_n))
  
  kept <- x %>%
    filter(n >= thr) %>%
    arrange(desc(n), desc(nchar(term)), term) %>%
    pull(term)
  
  if (length(kept) == 0) {
    # fallback: take the top term
    kept <- x %>% arrange(desc(n), desc(nchar(term)), term) %>% slice_head(n = 1) %>% pull(term)
  }
  
  paste(unique(kept), collapse = collapse)
}

# condense the bag of words into a single annotation label, like a "key words" label
df_out2 <- df_out %>%
  mutate(single_annotation = purrr::map_chr(bow_terms, make_label_maxanchored, p = 0.5, min_n = 1))


#print heads of the dataframes for debugging
print("patents_filtered head:")
print(head(patents_filtered))
print("df_out2 head:")
print(head(df_out2))


out_df=non_patent_db_unique %>% 
  mutate(single_annotation = antigen) %>% 
  bind_rows(df_out2 %>% mutate( dataset='patents' ))   
  
#score lexical similarity between all the values in the single_annotation column


write_parquet(out_df, output_path)
