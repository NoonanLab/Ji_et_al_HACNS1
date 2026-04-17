.libPaths(c("/gpfs/gibbs/project/noonan/yj345/R/4.2", .libPaths()))

library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(clustree)


FL <- list(list("CFLE9_HFLE9", "CFLE10_HFLE10", "CFLE11_HFLE11", "CFLE12_HFLE12"))
HL <- list(list("CHLE10_HHLE10", "CHLE11_HHLE11", "CHLE12_HHLE12"))
PA <- list(list("CPA1E9_HPA1E9", "CPA1E10_HPA1E10", "CPA1E11_HPA1E11", "CPA2E9_HPA2E9", "CPA2E10_HPA2E10", "CPA2E11_HPA2E11", "CPAE12_HPAE12"))

base_path <- "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/"
saved_path <- "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/"

transcription_factors <- readLines("~/project/allTFs_mm.txt")

process_FL <- function(samples, base_path, saved_path) {
  for (pair in samples) {
    sample1 <- pair[[1]]
    sample2 <- pair[[2]]
    sample3 <- pair[[3]]
    sample4 <- pair[[4]]
    
    
    sample1_obj <- readRDS(paste0(base_path, sample1, "_integrated.rds"))
    sample2_obj <- readRDS(paste0(base_path, sample2, "_integrated.rds"))
    sample3_obj <- readRDS(paste0(base_path, sample3, "_integrated.rds"))
    sample4_obj <- readRDS(paste0(base_path, sample4, "_integrated.rds"))
    
    
    sample1_obj$tissue <- "FL"
    sample2_obj$tissue <- "FL"
    sample3_obj$tissue <- "FL"
    sample4_obj$tissue <- "FL"
    
    sample1_obj$stage <- "E9.5"
    sample2_obj$stage <- "E10.5"
    sample3_obj$stage <- "E11.5"
    sample4_obj$stage <- "E12.5"
    
    
    merged_obj <- merge(x = sample1_obj, y = c(sample2_obj, sample3_obj, sample4_obj))
    all.genes <- rownames(merged_obj)
    merged_obj <- SCTransform(merged_obj)
    merged_obj <- RunPCA(merged_obj, features = transcription_factors)
    merged_obj <- RunUMAP(merged_obj, dims = 1:50)
    
    integrated_obj <- RunHarmony(object = merged_obj, group.by.vars = 'stage', project.dim = FALSE) %>%
      SCTransform() %>% RunPCA(verbose = FALSE, features = transcription_factors) %>% RunUMAP(reduction = "harmony", dims = 1:50, verbose = FALSE) %>%
      FindNeighbors(dims = 1:50, verbose = FALSE) %>% FindClusters(verbose = FALSE)
    
    DefaultAssay(integrated_obj) <- "RNA"
    integrated_obj <- NormalizeData(integrated_obj)
    integrated_obj <- ScaleData(integrated_obj, features = all.genes)
  
    saveRDS(integrated_obj, paste0(saved_path, "FL_integrated.rds"))
  }
}

process_FL(FL, base_path, saved_path)


process_HL <- function(samples, base_path, saved_path) {
  for (pair in samples) {
    sample1 <- pair[[1]]
    sample2 <- pair[[2]]
    sample3 <- pair[[3]]
    
    
    sample1_obj <- readRDS(paste0(base_path, sample1, "_integrated.rds"))
    sample2_obj <- readRDS(paste0(base_path, sample2, "_integrated.rds"))
    sample3_obj <- readRDS(paste0(base_path, sample3, "_integrated.rds"))
   
    
    sample1_obj$tissue <- "HL"
    sample2_obj$tissue <- "HL"
    sample3_obj$tissue <- "HL"
    
    sample1_obj$stage <- "E10.5"
    sample2_obj$stage <- "E11.5"
    sample3_obj$stage <- "E12.5"
    
    
    merged_obj <- merge(x = sample1_obj, y = c(sample2_obj, sample3_obj))
    all.genes <- rownames(merged_obj)
    merged_obj <- SCTransform(merged_obj)
    merged_obj <- RunPCA(merged_obj, features = transcription_factors)
    merged_obj <- RunUMAP(merged_obj, dims = 1:50)
    
    integrated_obj <- RunHarmony(object = merged_obj, group.by.vars = 'stage', project.dim = FALSE) %>%
      SCTransform() %>% RunPCA(verbose = FALSE, features = transcription_factors) %>% RunUMAP(reduction = "harmony", dims = 1:50, verbose = FALSE) %>%
      FindNeighbors(dims = 1:50, verbose = FALSE) %>% FindClusters(verbose = FALSE)
   
    DefaultAssay(integrated_obj) <- "RNA"
    integrated_obj <- NormalizeData(integrated_obj)
    integrated_obj <- ScaleData(integrated_obj, features = all.genes)
    
    saveRDS(integrated_obj, paste0(saved_path, "HL_integrated.rds"))
  }
}

process_HL(HL, base_path, saved_path)

process_PA <- function(samples, base_path, saved_path) {
  for (pair in samples) {
    sample1 <- pair[[1]]
    sample2 <- pair[[2]]
    sample3 <- pair[[3]]
    sample4 <- pair[[4]]
    sample5 <- pair[[5]]
    sample6 <- pair[[6]]
    sample7 <- pair[[7]]
    
    sample1_obj <- readRDS(paste0(base_path, sample1, "_integrated.rds"))
    sample2_obj <- readRDS(paste0(base_path, sample2, "_integrated.rds"))
    sample3_obj <- readRDS(paste0(base_path, sample3, "_integrated.rds"))
    sample4_obj <- readRDS(paste0(base_path, sample4, "_integrated.rds"))
    sample5_obj <- readRDS(paste0(base_path, sample5, "_integrated.rds"))
    sample6_obj <- readRDS(paste0(base_path, sample6, "_integrated.rds"))
    sample7_obj <- readRDS(paste0(base_path, sample7, "_integrated.rds"))
    
    sample1_obj$tissue <- "PA1"
    sample2_obj$tissue <- "PA1"
    sample3_obj$tissue <- "PA1"
    sample4_obj$tissue <- "PA2"
    sample5_obj$tissue <- "PA2"
    sample6_obj$tissue <- "PA2"
    sample7_obj$tissue <- "PA"
    sample1_obj$stage <- "E9.5"
    sample2_obj$stage <- "E10.5"
    sample3_obj$stage <- "E11.5"
    sample4_obj$stage <- "E9.5"
    sample5_obj$stage <- "E10.5"
    sample6_obj$stage <- "E11.5"
    sample7_obj$stage <- "E12.5"
    
    merged_obj <- merge(x = sample1_obj, y = c(sample2_obj, sample3_obj, sample4_obj, sample5_obj, sample6_obj, sample7_obj))
    all.genes <- rownames(merged_obj)
    merged_obj <- SCTransform(merged_obj)
    merged_obj <- RunPCA(merged_obj, features = transcription_factors)
    merged_obj <- RunUMAP(merged_obj, dims = 1:50)
    
    integrated_obj <- RunHarmony(object = merged_obj, group.by.vars = 'stage', project.dim = FALSE) %>%
      SCTransform() %>% RunPCA(verbose = FALSE, features = transcription_factors) %>% RunUMAP(reduction = "harmony", dims = 1:50, verbose = FALSE) %>%
      FindNeighbors(dims = 1:50, verbose = FALSE) %>% FindClusters(verbose = FALSE)
   
    DefaultAssay(integrated_obj) <- "RNA"
    integrated_obj <- NormalizeData(integrated_obj)
    integrated_obj <- ScaleData(integrated_obj, features = all.genes)
    
    saveRDS(integrated_obj, paste0(saved_path, "PA_integrated.rds"))
  }
}

process_PA(PA, base_path, saved_path)

rds_files <- list.files(path = saved_path, pattern = "\\.rds$", full.names = TRUE)

for (file in rds_files) {
  
  seurat_obj <- readRDS(file)
  
  
  resolutions <- seq(0, 1, by = 0.1)
  for (res in resolutions) {
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
  }
  
  
  clustree_plot <- clustree(seurat_obj, prefix = "SCT_snn_res.")
  
  
  plot_filename <- paste0(tools::file_path_sans_ext(basename(file)), "_clustree.png")
  plot_filepath <- file.path(saved_path, plot_filename)
  
  
  ggsave(filename = plot_filepath, plot = clustree_plot, width = 10, height = 10)
  
}
