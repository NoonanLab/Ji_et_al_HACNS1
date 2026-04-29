library(readr)
library(dplyr)
library(tidyr)
library(stringr)

GRN_BASE <- "/vast/palmer/scratch/noonan/yj345/snRNA/grn"
TOBIAS_BASE <- "/vast/palmer/scratch/noonan/yj345/ATAC/TOBIAS"

tissues <- c("FL", "HL", "PA1", "PA2")

tobias_tissue_map <- c(
  FL  = "FL",
  HL  = "HL",
  PA1 = "P1",
  PA2 = "P2"
)

filter_direct_targets <- function(tissue) {
  tobias_tissue <- tobias_tissue_map[[tissue]]
  
  grn_file <- file.path(
    GRN_BASE,
    tissue,
    "adj_sub",
    paste0("adj_", tissue, "_consensus_40runs_95pct.tsv")
  )
  
  adj_file <- file.path(
    TOBIAS_BASE,
    tobias_tissue,
    "adjacency_mousecase.txt"
  )
  
  out_file <- file.path(
    GRN_BASE,
    tissue,
    "adj_magic_direct_targets_final.tsv"
  )
  
  grn <- read_tsv(grn_file, show_col_types = FALSE)
  
  adj <- read_tsv(adj_file, show_col_types = FALSE)
  
  adj_long <- adj %>%
    separate_rows(Targets, sep = ",\\s*") %>%
    rename(
      TF = Source,
      target = Targets
    ) %>%
    mutate(
      TF = str_trim(TF),
      target = str_trim(target)
    ) %>%
    distinct(TF, target)
  
  direct_grn <- grn %>%
    mutate(
      TF = str_trim(TF),
      target = str_trim(target)
    ) %>%
    semi_join(adj_long, by = c("TF", "target"))
  
  write_tsv(direct_grn, out_file)
  
  message("Saved: ", out_file)
}

for (tissue in tissues) {
  filter_direct_targets(tissue)
}
