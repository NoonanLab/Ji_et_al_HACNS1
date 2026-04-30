library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(dplyr)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

base_dir <- "/vast/palmer/scratch/noonan/yj345/ATAC/new_E10"

tissues <- c("P1", "P2", "FL", "HL")

for (tissue in tissues) {
  
  peak_file <- file.path(base_dir, paste0(tissue, "_E10_Merged_peaks.bed"))
  out_file  <- file.path(base_dir, paste0(tissue, "_E10_Merged_peaks_annotated.bed"))
  
  peak_anno <- annotatePeak(
    peak_file,
    TxDb = txdb,
    tssRegion = c(-2000, 500),
    annoDb = "org.Mm.eg.db"
  )
  
  anno_df <- as.data.frame(peak_anno)
  
  annotated_bed <- anno_df %>%
    mutate(
      promoter_target_gene = ifelse(
        grepl("^Promoter", annotation),
        SYMBOL,
        NA
      )
    ) %>%
    select(
      seqnames,
      start,
      end,
      promoter_target_gene
    )
  
  write.table(
  annotated_bed,
  out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
  )
}
