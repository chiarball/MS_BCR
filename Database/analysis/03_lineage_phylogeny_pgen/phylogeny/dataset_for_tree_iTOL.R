
library(data.table)
library(ape, lib.loc = "/doctorai/chiarba/lib")

cluster_clone <- fread("/doctorai/chiarba/analysis/hilary/output/ALL_rows_with_clone_and_cluster_h1.tsv")

setDT(cluster_clone)   

## 1) Function to build iTOL datasets for one clone


make_itol_datasets_for_clone <- function(clone_of_interest,
                                         tree_dir = "/doctorai/chiarba/analysis/hilary/raxml",
                                         out_dir  = "/doctorai/chiarba/analysis/hilary/raxml/metadata_itol",
                                         seq_col  = "cdr3_aa",
                                         delta_outer = 8,
                                         delta_inner = 4,
                                         min_size = 10,
                                         max_size = 30) {
  
  cat("\n=== Processing clone", clone_of_interest, "===\n")
  
  ## 1) Read tree
  tree_file <- file.path(tree_dir, paste0("RAxML_bestTree.output_", clone_of_interest))
  if (!file.exists(tree_file)) {
    warning("Tree not found for clone ", clone_of_interest, ": ", tree_file)
    return(invisible(NULL))
  }
  
  tree <- read.tree(tree_file)
  tips <- data.table(node_id = tree$tip.label)
  tips[, sequence_id := gsub(" ", "_", node_id)]
  
  ## 2) Filter rows for this clone
  cc <- cluster_clone[clone_id == clone_of_interest]
  if (nrow(cc) == 0) {
    warning("Clone not found in cluster_clone: ", clone_of_interest)
    return(invisible(NULL))
  }
  
  ## count nt occurrences per sequence
  freq_nt <- cc[, .(n_nt = .N), by = sequence]
  
  ## tissues per sequence
  seq_tissue <- cc[, .(
    tissues   = paste(sort(unique(tissue)), collapse = ";"),
    n_tissues = uniqueN(tissue)
  ), by = sequence]
  
  ## basic metadata for input to iTOL
  meta_seq <- unique(cc[, .(
    sequence_id,
    sequence,
    seq_val = get(seq_col)
  )])
  
  meta_seq <- merge(meta_seq, freq_nt, by = "sequence", all.x = TRUE)
  meta_seq[is.na(n_nt), n_nt := 1L]
  
  meta_seq <- merge(meta_seq, seq_tissue, by = "sequence", all.x = TRUE)
  
  ## 3) Join with tree tips
  annot <- merge(tips, meta_seq, by = "sequence_id", all.x = FALSE)
  
  ## inner size based on n_nt
  if (length(unique(annot$n_nt)) == 1) {
    annot[, size_inner := (min_size + max_size) / 2]
  } else {
    annot[, size_inner := min_size +
             (n_nt - min(n_nt)) /
             (max(n_nt) - min(n_nt)) *
             (max_size - min_size)]
  }
  annot[, size_inner := round(size_inner, 1)]
  
  ## color palette for seq_val
  uniq_seq <- sort(unique(annot$seq_val))
  palette_fun <- colorRampPalette(c(
    "#007AFF", "#FF3B30", "#34C759", "#FF9500",
    "#AF52DE", "#FF2D55", "#5AC8FA", "#FFD60A", "#00C7BE"
  ))
  base_cols <- palette_fun(length(uniq_seq))
  names(base_cols) <- uniq_seq
  
  annot[, color_inner := paste0(base_cols[seq_val], "80")]
  
  ## DATASET 1: inner colored circles
 
  out_inner <- file.path(out_dir, paste0("clone", clone_of_interest, "_SYMBOL_bySequence.txt"))
  con1 <- file(out_inner, "w")
  
  writeLines("DATASET_SYMBOL", con1)
  writeLines("SEPARATOR TAB", con1)
  writeLines("", con1)
  writeLines(paste0("DATASET_LABEL\tClone", clone_of_interest, "_bySequence"), con1)
  writeLines("COLOR\t#000000", con1)
  writeLines(paste0("MAXIMUM_SIZE\t", max(annot$size_inner) + 12), con1)
  writeLines("", con1)
  writeLines("DATA", con1)
  
  lines_inner <- paste(
    annot$node_id,
    2,                      # circle
    annot$size_inner,
    annot$color_inner,
    1,                      # filled
    1.0,
    "",
    sep = "\t"
  )
  
  writeLines(lines_inner, con1)
  close(con1)
  cat("  -> inner symbols:", out_inner, "\n")
  
 
  ## shared nodes (>= 2 tissues)

  annot_shared <- annot[n_tissues >= 2]
  cat("  shared tips (>=2 tissues):", nrow(annot_shared), "\n")
  if (nrow(annot_shared) == 0) return(invisible(NULL))
  
  ## DATASET 2: outer red disk (background of ring)
 
  annot_shared[, size_outer := size_inner + delta_outer]
  
  out_outer <- file.path(out_dir, paste0("clone", clone_of_interest, "_SYMBOL_shared_border_outer.txt"))
  con_outer <- file(out_outer, "w")
  
  writeLines("DATASET_SYMBOL", con_outer)
  writeLines("SEPARATOR TAB", con_outer)
  writeLines("", con_outer)
  writeLines(paste0("DATASET_LABEL\tClone", clone_of_interest, "_shared_border_outer"), con_outer)
  writeLines("COLOR\t#FF0000", con_outer)
  writeLines(paste0("MAXIMUM_SIZE\t", max(annot_shared$size_outer) + 5), con_outer)
  writeLines("", con_outer)
  writeLines("DATA", con_outer)
  
  lines_outer <- paste(
    annot_shared$node_id,
    2,
    annot_shared$size_outer,
    "#FF0000FF",
    1,                      # filled
    1.0,
    annot_shared$tissues,
    sep = "\t"
  )
  
  writeLines(lines_outer, con_outer)
  close(con_outer)
  cat("  -> outer red ring base:", out_outer, "\n")
  
  ## DATASET 3: inner white mask
 
  annot_shared[, size_white := size_inner + delta_inner]
  
  out_white <- file.path(out_dir, paste0("clone", clone_of_interest, "_SYMBOL_shared_border_white.txt"))
  con_white <- file(out_white, "w")
  
  writeLines("DATASET_SYMBOL", con_white)
  writeLines("SEPARATOR TAB", con_white)
  writeLines("", con_white)
  writeLines(paste0("DATASET_LABEL\tClone", clone_of_interest, "_shared_border_white"), con_white)
  writeLines("COLOR\t#FFFFFF", con_white)
  writeLines(paste0("MAXIMUM_SIZE\t", max(annot_shared$size_white) + 5), con_white)
  writeLines("", con_white)
  writeLines("DATA", con_white)
  
  lines_white <- paste(
    annot_shared$node_id,
    2,
    annot_shared$size_white,
    "#FFFFFF",
    1,                      # filled
    1.0,
    "",
    sep = "\t"
  )
  
  writeLines(lines_white, con_white)
  close(con_white)
  cat("  -> inner white mask:", out_white, "\n")
  
  cat("=== DONE:", clone_of_interest, "===\n")
  invisible(NULL)
}

## 2) Run it for all desired clones

# list of clone_ids you want trees + iTOL files for:
clones_to_process <- c(40903, 19067)

lapply(clones_to_process, make_itol_datasets_for_clone)
