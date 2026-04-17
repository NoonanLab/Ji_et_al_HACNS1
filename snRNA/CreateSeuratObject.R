library(Seurat)
library(Matrix)

# Function to process a single sample
process_sample <- function(sample_path, project_name, output_dir) {
  data_dir <- file.path(sample_path, "Solo.out", "GeneFull", "filtered")
  sample_data <- Read10X(data.dir = data_dir)
  
  sample_name <- basename(sample_path)
  seurat_obj <- CreateSeuratObject(counts = sample_data, project = project_name, min.cells = 3, min.features = 200)
  
  spliced_counts <- readMM(file.path(sample_path, "Solo.out", "Velocyto", "filtered", "spliced.mtx"))
  unspliced_counts <- readMM(file.path(sample_path, "Solo.out", "Velocyto", "filtered", "unspliced.mtx"))
  barcodes <- readLines(file.path(sample_path, "Solo.out", "Velocyto", "filtered", "barcodes.tsv"))
  features <- readLines(file.path(sample_path, "Solo.out", "Velocyto", "filtered", "features.tsv"))
  
  rownames(spliced_counts) <- features
  rownames(unspliced_counts) <- features
  colnames(spliced_counts) <- barcodes
  colnames(unspliced_counts) <- barcodes
  
  cell_barcodes <- colnames(seurat_obj)
  common_barcodes <- intersect(colnames(spliced_counts), cell_barcodes)
  
  spliced_counts_subset <- spliced_counts[, common_barcodes]
  unspliced_counts_subset <- unspliced_counts[, common_barcodes]
  
  seurat_obj[["spliced"]] <- CreateAssayObject(counts = spliced_counts_subset)
  seurat_obj[["unspliced"]] <- CreateAssayObject(counts = unspliced_counts_subset)
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  saveRDS(seurat_obj, file.path(output_dir, paste0(sample_name, ".rds")))
}

# Base directory for all samples
base_dir <- "/vast/palmer/scratch/noonan/yj345/pipVelocyto"
output_dir <- "/vast/palmer/scratch/noonan/yj345/new"

# List of all samples
sample_paths <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
sample_paths <- sample_paths[grep("_Velocyto$", sample_paths)]

# Process each sample
for (sample_path in sample_paths) {
  process_sample(sample_path, project_name = "PA", output_dir = output_dir)
}
