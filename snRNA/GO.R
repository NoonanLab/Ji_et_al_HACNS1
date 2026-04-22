library(org.Mm.eg.db)
library(enrichplot)
library(GOSemSim)
library(dplyr)
library(readr)

ego <- enrichGO(
  gene          = genes,
  universe      = bg,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "ALL",       
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

df <- as.data.frame(ego)
