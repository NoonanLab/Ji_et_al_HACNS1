suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(URD))
library(Seurat)

FL <- readRDS("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_mes.rds")
HL <- readRDS("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_mes.rds")
PA_integrated <- readRDS("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA_integrated.rds")

PA1 <- PA_integrated[ ,PA_integrated$tissue %in% c("PA1", "PA")]
PA2 <- PA_integrated[ ,PA_integrated$tissue %in% "PA2"]

knitr::opts_chunk$set(echo = TRUE)
rgl::setupKnitr()

DefaultAssay(FL) <- "RNA"
FL_count <- GetAssayData(FL, assay = "RNA", layer = "counts")
FL_meta <- FL@meta.data
FL_URD <- createURD(count.data = FL_count, meta = FL_meta, min.cells = 3, min.counts = 3)

DefaultAssay(HL) <- "RNA"
HL_count <- GetAssayData(HL, assay = "RNA", layer = "counts")
HL_meta <- HL@meta.data
HL_URD <- createURD(count.data = HL_count, meta = HL_meta, min.cells = 3, min.counts = 3)

DefaultAssay(PA1) <- "RNA"
PA1_count <- GetAssayData(PA1, assay = "RNA", layer = "counts")
PA1_meta <- PA1@meta.data
PA1_URD <- createURD(count.data = PA1_count, meta = PA1_meta, min.cells = 3, min.counts = 3)

DefaultAssay(PA2) <- "RNA"
PA2_count <- GetAssayData(PA2, assay = "RNA", layer = "counts")
PA2_meta <- PA2@meta.data
PA2_URD <- createURD(count.data = PA2_count, meta = PA2_meta, min.cells = 3, min.counts = 3)

FL_URD@group.ids$stage <- as.character(FL_URD@meta[rownames(FL_URD@group.ids), "stage"])
HL_URD@group.ids$stage <- as.character(HL_URD@meta[rownames(HL_URD@group.ids), "stage"])
PA1_URD@group.ids$stage <- as.character(PA1_URD@meta[rownames(PA1_URD@group.ids), "stage"])
PA2_URD@group.ids$stage <- as.character(PA2_URD@meta[rownames(PA2_URD@group.ids), "stage"])

FL_stages <- sort(unique(FL_URD@group.ids$stage))
FL_var.by.stage <- lapply(1:(length(FL_stages) - 1), function(n) {
    findVariableGenes(
        FL_URD,
        cells.fit = cellsInCluster(FL_URD, "stage", FL_stages[n:(n + 1)]),
        set.object.var.genes = FALSE,
        diffCV.cutoff = 0.3,
        mean.min = 0.005,
        mean.max = 100,
        main.use = paste0("Stages ", FL_stages[n], " to ", FL_stages[n + 1]),
        do.plot = TRUE
    )
})
FL_var.genes <- sort(unique(unlist(FL_var.by.stage)))
FL_URD@var.genes <- FL_var.genes

HL_stages <- sort(unique(HL_URD@group.ids$stage))
HL_var.by.stage <- lapply(1:(length(HL_stages) - 1), function(n) {
    findVariableGenes(
        HL_URD,
        cells.fit = cellsInCluster(HL_URD, "stage", HL_stages[n:(n + 1)]),
        set.object.var.genes = FALSE,
        diffCV.cutoff = 0.3,
        mean.min = 0.005,
        mean.max = 100,
        main.use = paste0("Stages ", HL_stages[n], " to ", HL_stages[n + 1]),
        do.plot = TRUE
    )
})
HL_var.genes <- sort(unique(unlist(HL_var.by.stage)))
HL_URD@var.genes <- HL_var.genes

PA1_stages <- sort(unique(PA1_URD@group.ids$stage))
PA1_var.by.stage <- lapply(1:(length(PA1_stages) - 1), function(n) {
    findVariableGenes(
        PA1_URD,
        cells.fit = cellsInCluster(PA1_URD, "stage", PA1_stages[n:(n + 1)]),
        set.object.var.genes = FALSE,
        diffCV.cutoff = 0.3,
        mean.min = 0.005,
        mean.max = 100,
        main.use = paste0("Stages ", PA1_stages[n], " to ", PA1_stages[n + 1]),
        do.plot = TRUE
    )
})
PA1_var.genes <- sort(unique(unlist(PA1_var.by.stage)))
PA1_URD@var.genes <- PA1_var.genes

PA2_stages <- sort(unique(PA2_URD@group.ids$stage))
PA2_var.by.stage <- lapply(1:(length(PA2_stages) - 1), function(n) {
    findVariableGenes(
        PA2_URD,
        cells.fit = cellsInCluster(PA2_URD, "stage", PA2_stages[n:(n + 1)]),
        set.object.var.genes = FALSE,
        diffCV.cutoff = 0.3,
        mean.min = 0.005,
        mean.max = 100,
        main.use = paste0("Stages ", PA2_stages[n], " to ", PA2_stages[n + 1]),
        do.plot = TRUE
    )
})
PA2_var.genes <- sort(unique(unlist(PA2_var.by.stage)))
PA2_URD@var.genes <- PA2_var.genes

write(FL_var.genes,  file = "palmer_scratch/URD/FL/var_genes/var_genes.txt")
write(HL_var.genes,  file = "palmer_scratch/URD/HL/var_genes/var_genes.txt")
write(PA1_var.genes, file = "palmer_scratch/URD/PA1/var_genes/var_genes.txt")
write(PA2_var.genes, file = "palmer_scratch/URD/PA2/var_genes/var_genes.txt")

FL_URD <- calcPCA(FL_URD, mp.factor = 2)
pcSDPlot(FL_URD)
set.seed(19)
FL_URD <- calcTsne(FL_URD)
FL_URD <- graphClustering(FL_URD, dim.use = "pca", num.nn = c(15, 20, 30), do.jaccard = TRUE, method = "Louvain")
plotDim(FL_URD, "stage", plot.title = "FL stage")
plotDim(FL_URD, "cellcluster", plot.title = "FL cell type")

HL_URD <- calcPCA(HL_URD, mp.factor = 2)
pcSDPlot(HL_URD)
set.seed(19)
HL_URD <- calcTsne(HL_URD)
HL_URD <- graphClustering(HL_URD, dim.use = "pca", num.nn = c(15, 20, 30), do.jaccard = TRUE, method = "Louvain")
plotDim(HL_URD, "stage", plot.title = "HL stage")
plotDim(HL_URD, "cellcluster", plot.title = "HL cell type")

PA1_URD <- calcPCA(PA1_URD, mp.factor = 2)
pcSDPlot(PA1_URD)
set.seed(19)
PA1_URD <- calcTsne(PA1_URD)
PA1_URD <- graphClustering(PA1_URD, dim.use = "pca", num.nn = c(15, 20, 30), do.jaccard = TRUE, method = "Louvain")
plotDim(PA1_URD, "stage", plot.title = "PA1 stage")
plotDim(PA1_URD, "cellcluster", plot.title = "PA1 cell type")

PA2_URD <- calcPCA(PA2_URD, mp.factor = 2)
pcSDPlot(PA2_URD)
set.seed(19)
PA2_URD <- calcTsne(PA2_URD)
PA2_URD <- graphClustering(PA2_URD, dim.use = "pca", num.nn = c(15, 20, 30), do.jaccard = TRUE, method = "Louvain")
plotDim(PA2_URD, "stage", plot.title = "PA2 stage")
plotDim(PA2_URD, "cellcluster", plot.title = "PA2 cell type")

FL_URD <- calcKNN(FL_URD, nn = 100)
FL_outliers <- knnOutliers(FL_URD, nn.1 = 1, nn.2 = 20, x.max = 40,
                           slope.r = 1.1, int.r = 2.9, slope.b = 0.85, int.b = 10,
                           title = "FL Outliers")
FL_URD <- calcDM(FL_URD, knn = 170, sigma.use = 16)
plotDimArray(FL_URD, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "FL Diffusion Map", label = "stage", plot.title = "", legend = FALSE)
plotDim(FL_URD, "stage", reduction.use = "dm", transitions.plot = 10000, plot.title = "FL stage")

HL_URD <- calcKNN(HL_URD, nn = 100)
HL_outliers <- knnOutliers(HL_URD, nn.1 = 1, nn.2 = 20, x.max = 40,
                           slope.r = 1.1, int.r = 2.9, slope.b = 0.85, int.b = 10,
                           title = "HL Outliers")
HL_URD <- calcDM(HL_URD, knn = 170, sigma.use = 16)
plotDimArray(HL_URD, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "HL Diffusion Map", label = "stage", plot.title = "", legend = FALSE)
plotDim(HL_URD, "stage", reduction.use = "dm", transitions.plot = 10000, plot.title = "HL stage")

PA1_URD <- calcKNN(PA1_URD, nn = 100)
PA1_outliers <- knnOutliers(PA1_URD, nn.1 = 1, nn.2 = 20, x.max = 40,
                            slope.r = 1.1, int.r = 2.9, slope.b = 0.85, int.b = 10,
                            title = "PA1 Outliers")
PA1_URD <- calcDM(PA1_URD, knn = 170, sigma.use = 16)
plotDimArray(PA1_URD, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "PA1 Diffusion Map", label = "stage", plot.title = "", legend = FALSE)
plotDim(PA1_URD, "stage", reduction.use = "dm", transitions.plot = 10000, plot.title = "PA1 stage")

PA2_URD <- calcKNN(PA2_URD, nn = 100)
PA2_outliers <- knnOutliers(PA2_URD, nn.1 = 1, nn.2 = 20, x.max = 40,
                            slope.r = 1.1, int.r = 2.9, slope.b = 0.85, int.b = 10,
                            title = "PA2 Outliers")
PA2_URD <- calcDM(PA2_URD, knn = 170, sigma.use = 16)
plotDimArray(PA2_URD, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "PA2 Diffusion Map", label = "stage", plot.title = "", legend = FALSE)
plotDim(PA2_URD, "stage", reduction.use = "dm", transitions.plot = 10000, plot.title = "PA2 stage")

FL_root.cells <- rownames(FL_URD@meta)[FL_URD@meta$stage == "E9.5"]
HL_root.cells <- rownames(HL_URD@meta)[HL_URD@meta$stage == "E10.5"]
PA1_root.cells <- rownames(PA1_URD@meta)[PA1_URD@meta$stage == "E9.5"]
PA2_root.cells <- rownames(PA2_URD@meta)[PA2_URD@meta$stage == "E9.5"]

FL_floods <- floodPseudotime(FL_URD, root.cells = FL_root.cells, n = 150, minimum.cells.flooded = 2, verbose = FALSE)
FL_URD <- floodPseudotimeProcess(FL_URD, FL_floods, floods.name = "pseudotime", max.frac.NA = 0.4, pseudotime.fun = mean, stability.div = 20)

HL_floods <- floodPseudotime(HL_URD, root.cells = HL_root.cells, n = 150, minimum.cells.flooded = 2, verbose = FALSE)
HL_URD <- floodPseudotimeProcess(HL_URD, HL_floods, floods.name = "pseudotime", max.frac.NA = 0.4, pseudotime.fun = mean, stability.div = 20)

PA1_floods <- floodPseudotime(PA1_URD, root.cells = PA1_root.cells, n = 150, minimum.cells.flooded = 2, verbose = FALSE)
PA1_URD <- floodPseudotimeProcess(PA1_URD, PA1_floods, floods.name = "pseudotime", max.frac.NA = 0.4, pseudotime.fun = mean, stability.div = 20)

PA2_floods <- floodPseudotime(PA2_URD, root.cells = PA2_root.cells, n = 150, minimum.cells.flooded = 2, verbose = FALSE)
PA2_URD <- floodPseudotimeProcess(PA2_URD, PA2_floods, floods.name = "pseudotime", max.frac.NA = 0.4, pseudotime.fun = mean, stability.div = 20)

pseudotimePlotStabilityOverall(FL_URD)
plotDim(FL_URD, "pseudotime")
plotDists(FL_URD, "pseudotime", "stage", plot.title = "FL pseudotime by stage")

pseudotimePlotStabilityOverall(HL_URD)
plotDim(HL_URD, "pseudotime")
plotDists(HL_URD, "pseudotime", "stage", plot.title = "HL pseudotime by stage")

pseudotimePlotStabilityOverall(PA1_URD)
plotDim(PA1_URD, "pseudotime")
plotDists(PA1_URD, "pseudotime", "stage", plot.title = "PA1 pseudotime by stage")

pseudotimePlotStabilityOverall(PA2_URD)
plotDim(PA2_URD, "pseudotime")
plotDists(PA2_URD, "pseudotime", "stage", plot.title = "PA2 pseudotime by stage")

