library(Seurat)
library(ggplot2)
library(cellAlign)
library(reshape2)
library(pheatmap)
library(SeuratData)
library(SeuratDisk)

Convert(""/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_Cartilage.h5ad"", dest = "h5seurat", overwrite = TRUE)
FL_Cartilage <- LoadH5Seurat("FL_Cartilage.h5seurat")
FL_Cartilage

Convert(""/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_Cartilage.h5ad"", dest = "h5seurat", overwrite = TRUE)
HL_Cartilage <- LoadH5Seurat("HL_Cartilage.h5seurat")
HL_Cartilage

Convert(""/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA_Osteoblasts.h5ad"", dest = "h5seurat", overwrite = TRUE)
PA_Osteoblasts <- LoadH5Seurat("PA_Osteoblasts.h5seurat")
PA_Osteoblasts

cellalign_lineage <- function(seurat_object) {
    DefaultAssay(seurat_object) <- "RNA"
    
    H_lineage <- seurat_object[ ,seurat_object$genotype %in% "Human"]
    C_lineage <- seurat_object[ ,seurat_object$genotype %in% "Chimp"]
    
    # Find variable features
    H_lineage <- FindVariableFeatures(H_lineage, nfeatures = 3000)
    C_lineage <- FindVariableFeatures(C_lineage, nfeatures = 3000)
    
    # Identify shared HVFs
    shared_HVFs <- union(VariableFeatures(H_lineage), VariableFeatures(C_lineage))
    
    # Filter to shared HVFs
    H_lineage_filtered <- H_lineage[rownames(H_lineage) %in% shared_HVFs,]
    C_lineage_filtered <- C_lineage[rownames(C_lineage) %in% shared_HVFs,]
    
    H_lineage_filtered <- GetAssayData(H_lineage_filtered, assay = "RNA", layer = "counts")
    C_lineage_filtered <- GetAssayData(C_lineage_filtered, assay = "RNA", layer = "counts")
    
    # Filter to non-zero features
    shared_non_zero_HVFs <- intersect(rownames(H_lineage_filtered[rowSums(H_lineage_filtered) > 0,]),
                                      rownames(C_lineage_filtered[rowSums(C_lineage_filtered) > 0,]))
    
    H_lineage_filtered <- H_lineage_filtered[shared_non_zero_HVFs,]
    C_lineage_filtered <- C_lineage_filtered[shared_non_zero_HVFs,]
    
    # Sort rows
    H_lineage_filtered <- H_lineage_filtered[sort(rownames(H_lineage_filtered)),]
    C_lineage_filtered <- C_lineage_filtered[sort(rownames(C_lineage_filtered)),]
    
    # Further alignment
    H_lineage_filtered <- H_lineage_filtered[rownames(H_lineage_filtered) %in% rownames(C_lineage_filtered),]
    C_lineage_filtered <- C_lineage_filtered[rownames(C_lineage_filtered) %in% rownames(H_lineage_filtered),]
    
    C_lineage_filtered <- C_lineage_filtered[rownames(H_lineage_filtered),]
    
    # Rank pseudotime
    H_lineage$ranked_pseudotime <- rank(H_lineage$slingshot_pseudotime)
    C_lineage$ranked_pseudotime <- rank(C_lineage$slingshot_pseudotime)
    
    trajHumRTN <- H_lineage$ranked_pseudotime / max(H_lineage$ranked_pseudotime)
    trajChRTN <- C_lineage$ranked_pseudotime / max(C_lineage$ranked_pseudotime)
    
    names(trajHumRTN) <- colnames(H_lineage)
    names(trajChRTN) <- colnames(C_lineage)
    
    trajHumRTN <- sort(trajHumRTN)
    trajChRTN <- sort(trajChRTN)
    
    H_lineage_filtered <- H_lineage_filtered[, names(trajHumRTN)]
    C_lineage_filtered <- C_lineage_filtered[, names(trajChRTN)]
    
    trajHumRTN <- trajHumRTN[colnames(H_lineage_filtered)]
    trajChRTN <- trajChRTN[colnames(C_lineage_filtered)]
    
    # Interpolation
    interGlobalHumRTN <- cellAlign::interWeights(expDataBatch = H_lineage_filtered,
                                                 trajCond = trajHumRTN,
                                                 winSz = 0.1, numPts = 200)
    interGlobalChRTN <- cellAlign::interWeights(expDataBatch = C_lineage_filtered,
                                                trajCond = trajChRTN,
                                                winSz = 0.1, numPts = 200)
    
    interScaledGlobalHumRTN <- cellAlign::scaleInterpolate(interGlobalHumRTN)
    interScaledGlobalChRTN <- cellAlign::scaleInterpolate(interGlobalChRTN)
    
    # Calculate distance matrix and plot heatmap
    interRTNDistMat <- calcDistMat(interScaledGlobalHumRTN$scaledData, interScaledGlobalChRTN$scaledData, dist.method = 'Euclidean')
    pheatmap(interRTNDistMat, cluster_cols = FALSE, cluster_rows = FALSE, main = "Human vs. Chimp Cell Distances", show_rownames = FALSE, show_colnames = FALSE, display_numbers = FALSE)
    
    # Alignment and plot
    interRTNAlign <- globalAlign(interScaledGlobalHumRTN$scaledData, interScaledGlobalChRTN$scaledData,
                                 scores = list(query = interScaledGlobalHumRTN$traj,
                                               ref = interScaledGlobalChRTN$traj),
                                 sigCalc = FALSE, numPerm = 30)
    plotAlign(interRTNAlign)
    
    # Mapping and plot
    mappingRTNAlign <- mapRealDataGlobal(interRTNAlign, intTrajQuery = interScaledGlobalHumRTN$traj, realTrajQuery = trajHumRTN,
                                         intTrajRef = interScaledGlobalChRTN$traj, realTrajRef = trajChRTN)
    plotMapping(mappingRTNAlign)
    
}

cellalign_lineage(FL_Cartilage)
cellalign_lineage(HL_Cartilage)
cellalign_lineage(PA_Osteoblasts)
