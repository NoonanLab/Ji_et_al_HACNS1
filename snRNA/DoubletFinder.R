library(Seurat)
library(dplyr)
library(cowplot)
library(DoubletFinder)

out_dir <- "~/palmer_scratch/new"

samples <- c(
  "CPA1E9","CPA2E9","HPA1E9","HPA2E9","CFLE9","HFLE9",
  "CPA1E10","CPA2E10","HPA1E10","HPA2E10","CFLE10","CHLE10","HFLE10","HHLE10",
  "CPA1E11","CPA2E11","HPA1E11","HPA2E11","CFLE11","CHLE11","HFLE11","HHLE11",
  "CPAE12","HPAE12","CFLE12","CHLE12","HFLE12","HHLE12"
)

for (sample in samples) {
  obj <- get(paste0(sample, "_Velocyto_CCR"))
  assign(
    paste0("sweep.res.list_", sample),
    paramSweep_v3(obj, PCs = 1:30, sct = TRUE),
    envir = .GlobalEnv
  )
}

for (sample in samples) {
  sweep.res <- get(paste0("sweep.res.list_", sample))
  assign(
    paste0("sweep.stats_", sample),
    summarizeSweep(sweep.res, GT = FALSE),
    envir = .GlobalEnv
  )
}

for (sample in samples) {
  sweep.stats <- get(paste0("sweep.stats_", sample))
  assign(
    paste0("bcmvn_", sample),
    find.pK(sweep.stats),
    envir = .GlobalEnv
  )
}

for (sample in samples) {
  obj <- get(paste0(sample, "_Velocyto_CCR"))
  assign(
    paste0("annotations_", sample),
    obj@meta.data$seurat_clusters,
    envir = .GlobalEnv
  )
}

for (sample in samples) {
  annotations <- get(paste0("annotations_", sample))
  assign(
    paste0("homotypic.prop_", sample),
    modelHomotypic(annotations),
    envir = .GlobalEnv
  )
}

nExp_poi_CPA1E9 <- round(0.075 * nrow(CPA1E9_Velocyto_CCR@meta.data))
nExp_poi_CPA2E9 <- round(0.075 * nrow(CPA2E9_Velocyto_CCR@meta.data))
nExp_poi_HPA1E9 <- round(0.10  * nrow(HPA1E9_Velocyto_CCR@meta.data))
nExp_poi_HPA2E9 <- round(0.075 * nrow(HPA2E9_Velocyto_CCR@meta.data))
nExp_poi_CFLE9  <- round(0.075 * nrow(CFLE9_Velocyto_CCR@meta.data))
nExp_poi_HFLE9  <- round(0.075 * nrow(HFLE9_Velocyto_CCR@meta.data))

nExp_poi_CPA1E10 <- round(0.075 * nrow(CPA1E10_Velocyto_CCR@meta.data))
nExp_poi_CPA2E10 <- round(0.075 * nrow(CPA2E10_Velocyto_CCR@meta.data))
nExp_poi_HPA1E10 <- round(0.075 * nrow(HPA1E10_Velocyto_CCR@meta.data))
nExp_poi_HPA2E10 <- round(0.075 * nrow(HPA2E10_Velocyto_CCR@meta.data))
nExp_poi_CFLE10  <- round(0.075 * nrow(CFLE10_Velocyto_CCR@meta.data))
nExp_poi_CHLE10  <- round(0.075 * nrow(CHLE10_Velocyto_CCR@meta.data))
nExp_poi_HFLE10  <- round(0.075 * nrow(HFLE10_Velocyto_CCR@meta.data))
nExp_poi_HHLE10  <- round(0.075 * nrow(HHLE10_Velocyto_CCR@meta.data))

nExp_poi_CPA1E11 <- round(0.075 * nrow(CPA1E11_Velocyto_CCR@meta.data))
nExp_poi_CPA2E11 <- round(0.075 * nrow(CPA2E11_Velocyto_CCR@meta.data))
nExp_poi_HPA1E11 <- round(0.075 * nrow(HPA1E11_Velocyto_CCR@meta.data))
nExp_poi_HPA2E11 <- round(0.075 * nrow(HPA2E11_Velocyto_CCR@meta.data))
nExp_poi_CFLE11  <- round(0.075 * nrow(CFLE11_Velocyto_CCR@meta.data))
nExp_poi_CHLE11  <- round(0.15  * nrow(CHLE11_Velocyto_CCR@meta.data))
nExp_poi_HFLE11  <- round(0.075 * nrow(HFLE11_Velocyto_CCR@meta.data))
nExp_poi_HHLE11  <- round(0.075 * nrow(HHLE11_Velocyto_CCR@meta.data))

nExp_poi_CPAE12 <- round(0.075 * nrow(CPAE12_Velocyto_CCR@meta.data))
nExp_poi_HPAE12 <- round(0.075 * nrow(HPAE12_Velocyto_CCR@meta.data))
nExp_poi_CFLE12 <- round(0.075 * nrow(CFLE12_Velocyto_CCR@meta.data))
nExp_poi_CHLE12 <- round(0.075 * nrow(CHLE12_Velocyto_CCR@meta.data))
nExp_poi_HFLE12 <- round(0.10  * nrow(HFLE12_Velocyto_CCR@meta.data))
nExp_poi_HHLE12 <- round(0.15  * nrow(HHLE12_Velocyto_CCR@meta.data))

for (sample in samples) {
  nExp_poi <- get(paste0("nExp_poi_", sample))
  homotypic.prop <- get(paste0("homotypic.prop_", sample))
  assign(
    paste0("nExp_poi.adj_", sample),
    round(nExp_poi * (1 - homotypic.prop)),
    envir = .GlobalEnv
  )
}

CPA1E9_DF <- doubletFinder_v3(CPA1E9_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi_CPA1E9, reuse.pANN = FALSE, sct = TRUE)
CPA2E9_DF <- doubletFinder_v3(CPA2E9_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CPA2E9, reuse.pANN = FALSE, sct = TRUE)
HPA1E9_DF <- doubletFinder_v3(HPA1E9_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.16, nExp = nExp_poi_HPA1E9, reuse.pANN = FALSE, sct = TRUE)
HPA2E9_DF <- doubletFinder_v3(HPA2E9_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_HPA2E9, reuse.pANN = FALSE, sct = TRUE)

CFLE9_DF <- doubletFinder_v3(CFLE9_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CFLE9, reuse.pANN = FALSE, sct = TRUE)
HFLE9_DF <- doubletFinder_v3(HFLE9_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_HFLE9, reuse.pANN = FALSE, sct = TRUE)

CPA1E10_DF <- doubletFinder_v3(CPA1E10_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CPA1E10, reuse.pANN = FALSE, sct = TRUE)
CPA2E10_DF <- doubletFinder_v3(CPA2E10_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CPA2E10, reuse.pANN = FALSE, sct = TRUE)
HPA1E10_DF <- doubletFinder_v3(HPA1E10_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_HPA1E10, reuse.pANN = FALSE, sct = TRUE)
HPA2E10_DF <- doubletFinder_v3(HPA2E10_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.04, nExp = nExp_poi_HPA2E10, reuse.pANN = FALSE, sct = TRUE)

CFLE10_DF <- doubletFinder_v3(CFLE10_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CFLE10, reuse.pANN = FALSE, sct = TRUE)
CHLE10_DF <- doubletFinder_v3(CHLE10_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CHLE10, reuse.pANN = FALSE, sct = TRUE)
HFLE10_DF <- doubletFinder_v3(HFLE10_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.16, nExp = nExp_poi_HFLE10, reuse.pANN = FALSE, sct = TRUE)
HHLE10_DF <- doubletFinder_v3(HHLE10_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_HHLE10, reuse.pANN = FALSE, sct = TRUE)

CPA1E11_DF <- doubletFinder_v3(CPA1E11_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CPA1E11, reuse.pANN = FALSE, sct = TRUE)
CPA2E11_DF <- doubletFinder_v3(CPA2E11_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CPA2E11, reuse.pANN = FALSE, sct = TRUE)
HPA1E11_DF <- doubletFinder_v3(HPA1E11_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_HPA1E11, reuse.pANN = FALSE, sct = TRUE)
HPA2E11_DF <- doubletFinder_v3(HPA2E11_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_HPA2E11, reuse.pANN = FALSE, sct = TRUE)

CFLE11_DF <- doubletFinder_v3(CFLE11_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CFLE11, reuse.pANN = FALSE, sct = TRUE)
CHLE11_DF <- doubletFinder_v3(CHLE11_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.07, nExp = nExp_poi_CHLE11, reuse.pANN = FALSE, sct = TRUE)
HFLE11_DF <- doubletFinder_v3(HFLE11_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_HFLE11, reuse.pANN = FALSE, sct = TRUE)
HHLE11_DF <- doubletFinder_v3(HHLE11_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_HHLE11, reuse.pANN = FALSE, sct = TRUE)

CPAE12_DF <- doubletFinder_v3(CPAE12_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.04, nExp = nExp_poi_CPAE12, reuse.pANN = FALSE, sct = TRUE)
HPAE12_DF <- doubletFinder_v3(HPAE12_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.04, nExp = nExp_poi_HPAE12, reuse.pANN = FALSE, sct = TRUE)

CFLE12_DF <- doubletFinder_v3(CFLE12_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CFLE12, reuse.pANN = FALSE, sct = TRUE)
CHLE12_DF <- doubletFinder_v3(CHLE12_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.11, nExp = nExp_poi_CHLE12, reuse.pANN = FALSE, sct = TRUE)
HFLE12_DF <- doubletFinder_v3(HFLE12_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.21, nExp = nExp_poi_HFLE12, reuse.pANN = FALSE, sct = TRUE)
HHLE12_DF <- doubletFinder_v3(HHLE12_Velocyto_CCR, PCs = 1:30, pN = 0.25, pK = 0.29, nExp = nExp_poi_HHLE12, reuse.pANN = FALSE, sct = TRUE)

Idents(CPA1E9_DF) <- "DF.classifications_0.25_0.09_1874"
CPA1E9_Velocyto_CCR <- subset(CPA1E9_DF, ident = "Singlet")
saveRDS(CPA1E9_Velocyto_CCR, "~/palmer_scratch/new/CPA1E9_Velocyto_CCR.rds")

Idents(CPA2E9_DF) <- "DF.classifications_0.25_0.11_1900"
CPA2E9_Velocyto_CCR <- subset(CPA2E9_DF, ident = "Singlet")
saveRDS(CPA2E9_Velocyto_CCR, "~/palmer_scratch/new/CPA2E9_Velocyto_CCR.rds")

Idents(HPA1E9_DF) <- "DF.classifications_0.25_0.16_2500"
HPA1E9_Velocyto_CCR <- subset(HPA1E9_DF, ident = "Singlet")
saveRDS(HPA1E9_Velocyto_CCR, "~/palmer_scratch/new/HPA1E9_Velocyto_CCR.rds")

Idents(HPA2E9_DF) <- "DF.classifications_0.25_0.11_2100"
HPA2E9_Velocyto_CCR <- subset(HPA2E9_DF, ident = "Singlet")
saveRDS(HPA2E9_Velocyto_CCR, "~/palmer_scratch/new/HPA2E9_Velocyto_CCR.rds")

Idents(CFLE9_DF) <- "DF.classifications_0.25_0.11_1800"
CFLE9_Velocyto_CCR <- subset(CFLE9_DF, ident = "Singlet")
saveRDS(CFLE9_Velocyto_CCR, "~/palmer_scratch/new/CFLE9_Velocyto_CCR.rds")

Idents(HFLE9_DF) <- "DF.classifications_0.25_0.11_1850"
HFLE9_Velocyto_CCR <- subset(HFLE9_DF, ident = "Singlet")
saveRDS(HFLE9_Velocyto_CCR, "~/palmer_scratch/new/HFLE9_Velocyto_CCR.rds")

Idents(CPA1E10_DF) <- "DF.classifications_0.25_0.11_2000"
CPA1E10_Velocyto_CCR <- subset(CPA1E10_DF, ident = "Singlet")
saveRDS(CPA1E10_Velocyto_CCR, "~/palmer_scratch/new/CPA1E10_Velocyto_CCR.rds")

Idents(CPA2E10_DF) <- "DF.classifications_0.25_0.11_2050"
CPA2E10_Velocyto_CCR <- subset(CPA2E10_DF, ident = "Singlet")
saveRDS(CPA2E10_Velocyto_CCR, "~/palmer_scratch/new/CPA2E10_Velocyto_CCR.rds")

Idents(HPA1E10_DF) <- "DF.classifications_0.25_0.11_1980"
HPA1E10_Velocyto_CCR <- subset(HPA1E10_DF, ident = "Singlet")
saveRDS(HPA1E10_Velocyto_CCR, "~/palmer_scratch/new/HPA1E10_Velocyto_CCR.rds")

Idents(HPA2E10_DF) <- "DF.classifications_0.25_0.04_1760"
HPA2E10_Velocyto_CCR <- subset(HPA2E10_DF, ident = "Singlet")
saveRDS(HPA2E10_Velocyto_CCR, "~/palmer_scratch/new/HPA2E10_Velocyto_CCR.rds")

Idents(CFLE10_DF) <- "DF.classifications_0.25_0.11_2100"
CFLE10_Velocyto_CCR <- subset(CFLE10_DF, ident = "Singlet")
saveRDS(CFLE10_Velocyto_CCR, "~/palmer_scratch/new/CFLE10_Velocyto_CCR.rds")

Idents(CHLE10_DF) <- "DF.classifications_0.25_0.11_2200"
CHLE10_Velocyto_CCR <- subset(CHLE10_DF, ident = "Singlet")
saveRDS(CHLE10_Velocyto_CCR, "~/palmer_scratch/new/CHLE10_Velocyto_CCR.rds")

Idents(HFLE10_DF) <- "DF.classifications_0.25_0.16_2300"
HFLE10_Velocyto_CCR <- subset(HFLE10_DF, ident = "Singlet")
saveRDS(HFLE10_Velocyto_CCR, "~/palmer_scratch/new/HFLE10_Velocyto_CCR.rds")

Idents(HHLE10_DF) <- "DF.classifications_0.25_0.11_2150"
HHLE10_Velocyto_CCR <- subset(HHLE10_DF, ident = "Singlet")
saveRDS(HHLE10_Velocyto_CCR, "~/palmer_scratch/new/HHLE10_Velocyto_CCR.rds")

Idents(CPA1E11_DF) <- "DF.classifications_0.25_0.11_2400"
CPA1E11_Velocyto_CCR <- subset(CPA1E11_DF, ident = "Singlet")
saveRDS(CPA1E11_Velocyto_CCR, "~/palmer_scratch/new/CPA1E11_Velocyto_CCR.rds")

Idents(CPA2E11_DF) <- "DF.classifications_0.25_0.11_2350"
CPA2E11_Velocyto_CCR <- subset(CPA2E11_DF, ident = "Singlet")
saveRDS(CPA2E11_Velocyto_CCR, "~/palmer_scratch/new/CPA2E11_Velocyto_CCR.rds")

Idents(HPA1E11_DF) <- "DF.classifications_0.25_0.11_2450"
HPA1E11_Velocyto_CCR <- subset(HPA1E11_DF, ident = "Singlet")
saveRDS(HPA1E11_Velocyto_CCR, "~/palmer_scratch/new/HPA1E11_Velocyto_CCR.rds")

Idents(HPA2E11_DF) <- "DF.classifications_0.25_0.11_2380"
HPA2E11_Velocyto_CCR <- subset(HPA2E11_DF, ident = "Singlet")
saveRDS(HPA2E11_Velocyto_CCR, "~/palmer_scratch/new/HPA2E11_Velocyto_CCR.rds")

Idents(CFLE11_DF) <- "DF.classifications_0.25_0.11_2500"
CFLE11_Velocyto_CCR <- subset(CFLE11_DF, ident = "Singlet")
saveRDS(CFLE11_Velocyto_CCR, "~/palmer_scratch/new/CFLE11_Velocyto_CCR.rds")

Idents(CHLE11_DF) <- "DF.classifications_0.25_0.07_2777"
CHLE11_Velocyto_CCR <- subset(CHLE11_DF, ident = "Singlet")
saveRDS(CHLE11_Velocyto_CCR, "~/palmer_scratch/new/CHLE11_Velocyto_CCR.rds")

Idents(HFLE11_DF) <- "DF.classifications_0.25_0.11_2520"
HFLE11_Velocyto_CCR <- subset(HFLE11_DF, ident = "Singlet")
saveRDS(HFLE11_Velocyto_CCR, "~/palmer_scratch/new/HFLE11_Velocyto_CCR.rds")

Idents(HHLE11_DF) <- "DF.classifications_0.25_0.11_2480"
HHLE11_Velocyto_CCR <- subset(HHLE11_DF, ident = "Singlet")
saveRDS(HHLE11_Velocyto_CCR, "~/palmer_scratch/new/HHLE11_Velocyto_CCR.rds")

Idents(CPAE12_DF) <- "DF.classifications_0.25_0.04_1600"
CPAE12_Velocyto_CCR <- subset(CPAE12_DF, ident = "Singlet")
saveRDS(CPAE12_Velocyto_CCR, "~/palmer_scratch/new/CPAE12_Velocyto_CCR.rds")

Idents(HPAE12_DF) <- "DF.classifications_0.25_0.04_1650"
HPAE12_Velocyto_CCR <- subset(HPAE12_DF, ident = "Singlet")
saveRDS(HPAE12_Velocyto_CCR, "~/palmer_scratch/new/HPAE12_Velocyto_CCR.rds")

Idents(CFLE12_DF) <- "DF.classifications_0.25_0.11_1900"
CFLE12_Velocyto_CCR <- subset(CFLE12_DF, ident = "Singlet")
saveRDS(CFLE12_Velocyto_CCR, "~/palmer_scratch/new/CFLE12_Velocyto_CCR.rds")

Idents(CHLE12_DF) <- "DF.classifications_0.25_0.11_1950"
CHLE12_Velocyto_CCR <- subset(CHLE12_DF, ident = "Singlet")
saveRDS(CHLE12_Velocyto_CCR, "~/palmer_scratch/new/CHLE12_Velocyto_CCR.rds")

Idents(HFLE12_DF) <- "DF.classifications_0.25_0.21_2100"
HFLE12_Velocyto_CCR <- subset(HFLE12_DF, ident = "Singlet")
saveRDS(HFLE12_Velocyto_CCR, "~/palmer_scratch/new/HFLE12_Velocyto_CCR.rds")

Idents(HHLE12_DF) <- "DF.classifications_0.25_0.29_3200"
HHLE12_Velocyto_CCR <- subset(HHLE12_DF, ident = "Singlet")
saveRDS(HHLE12_Velocyto_CCR, "~/palmer_scratch/new/HHLE12_Velocyto_CCR.rds")
