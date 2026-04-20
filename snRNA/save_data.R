save_data_by_tissue <- function(seurat_obj, tissue, save_dir) {
  DefaultAssay(seurat_obj) <- "RNA"
  # Get RNA counts
  RNA_count <- GetAssayData(seurat_obj, slot = "counts")
  RNA_count <- as.matrix(RNA_count)
  RNA_count <- t(RNA_count)

  rna_count_file <- file.path(save_dir, paste0(tissue, "_RNA_counts.tsv"))
  metadata_file <- file.path(save_dir, paste0(tissue, "_metadata.csv"))
  harmony_file <- file.path(save_dir, paste0(tissue, "_harmony.csv"))
 
  write.table(RNA_count, file = rna_count_file, quote = FALSE, sep = '\t')
  write.csv(seurat_obj@meta.data, file = metadata_file)
  write.csv(seurat_obj@reductions$harmony@cell.embeddings, file = harmony_file)
}

save_data_by_tissue(FL_mes, "FL", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/Phate/FL/")
save_data_by_tissue(HL_mes, "FL", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/Phate/HL/")
save_data_by_tissue(PA_mes, "FL", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/Phate/PA/")
