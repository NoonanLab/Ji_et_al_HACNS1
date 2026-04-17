.libPaths(c("/gpfs/gibbs/project/noonan/yj345/R/4.2", .libPaths()))

setwd("/vast/palmer/scratch/noonan/yj345/new")

library(Seurat)
library(dplyr)
library(cowplot)
library(stringr)

s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)

rds_files <- list.files(pattern = "\\.rds$")

for (file in rds_files) {

    data <- readRDS(file)
    DefaultAssay(data) <- "RNA"
    data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    data <- RunPCA(data, features = c(s.genes, g2m.genes), reduction.name = "cc.pca")
    data$CC.Difference <- data$S.Score - data$G2M.Score
    data <- ScaleData(data, vars.to.regress = "CC.Difference", features = rownames(data))
    data <- RunPCA(data, features = c(s.genes, g2m.genes), reduction.name = "cc.pca")
    DefaultAssay(data) <- "SCT"
    new_filename <- paste0(tools::file_path_sans_ext(file), "_CCR.rds")
    saveRDS(data, new_filename)
    
}
