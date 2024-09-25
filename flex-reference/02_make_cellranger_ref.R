library(data.table)
library(dplyr)
library(Seurat)
library(readxl)
library(stringr)

# Oracle: https://kb.10xgenomics.com/hc/en-us/articles/17824038751885-How-to-Add-Custom-Probes-to-Cell-Ranger-and-Space-Ranger

# Import designed sequences
probe_ids_df <- read_excel("../HHV6_probes_with_human_background_1barcode_info.xlsx") 
probe_ids_df$gene <- stringr::str_split_fixed(probe_ids_df$name, pattern = "=|\\[|\\]", n = 6)[,3]
probe_ids_df$probe_u <- make.unique(probe_ids_df$gene)

n_probes <- dim(probe_ids_df)[1]
probe_ids_df$hash_id_ref <- paste0("ffff", str_pad(1:n_probes, 3, pad = "0"))

# create data frame

odf <- data.frame(gene = probe_ids_df$gene,
           probe_seq = paste0(probe_ids_df$lhs_seq, probe_ids_df$rhs_seq),
           probe_id = paste0(probe_ids_df$gene, "|", probe_ids_df$gene, "|", probe_ids_df$hash_id_ref),
           included = "TRUE", region = "unspliced")
odf
write.table(odf, file = "hhv6_custom_flex.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
