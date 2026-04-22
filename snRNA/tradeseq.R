library(Seurat)
library(tradeSeq)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cellAlign)
library(Matrix)
library(tibble)

FL <- readRDS("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_mes.rds")
HL <- readRDS("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_mes.rds")
PA_integrated <- readRDS("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA_integrated.rds")

run_tradeseq <- function(
  seu,
  pseudotime_col,
  lineage_name,
  outdir
) {
  hvf <- VariableFeatures(seu)
  
  H_HVFs <- seu[, seu$genotype %in% "Human"]
  C_HVFs <- seu[, seu$genotype %in% "Chimp"]
  
  H_HVFs <- FindVariableFeatures(H_HVFs, assay = "RNA", nfeatures = 3000)
  C_HVFs <- FindVariableFeatures(C_HVFs, assay = "RNA", nfeatures = 3000)
  
  shared_HVFs <- union(VariableFeatures(H_HVFs), VariableFeatures(C_HVFs))
  shared_HVFs <- union(shared_HVFs, hvf)
  
  all_genes <- shared_HVFs
  
  is_ribo <- grepl("^(Rps|Rpl|Mrps|Mrpl)", all_genes, perl = TRUE)
  is_hsp  <- grepl("^(Hsp|Hspa|Hspb|Hspc|Hspd|Hspe|Hsph|Dnaja|Dnajb|Dnajc)", all_genes, perl = TRUE)
  is_mt   <- grepl("^(MT-|mt-)", all_genes, perl = TRUE)
  
  shared_HVFs <- all_genes[!(is_ribo | is_hsp | is_mt)]
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
 
  lineage_subset <- seu[, !is.na(seu[[pseudotime_col, drop = TRUE]])]
  
  H_lineage <- lineage_subset[, lineage_subset$genotype %in% "Human"]
  C_lineage <- lineage_subset[, lineage_subset$genotype %in% "Chimp"]
  
  H_lineage_filtered <- GetAssayData(
    H_lineage[shared_HVFs, ],
    assay = "SCT",
    slot = "counts"
  )
  
  C_lineage_filtered <- GetAssayData(
    C_lineage[shared_HVFs, ],
    assay = "SCT",
    slot = "counts"
  )
  
  shared_non_zero_HVFs <- intersect(
    rownames(H_lineage_filtered[rowSums(H_lineage_filtered) > 0, ]),
    rownames(C_lineage_filtered[rowSums(C_lineage_filtered) > 0, ])
  )
  
  shared_non_zero_HVFs <- union(shared_non_zero_HVFs, "Gbx2")
  
  Lin_tradeSeq <- list()
  Lin_tradeSeq$counts <- as.matrix(
    lineage_subset@assays$RNA@counts[shared_non_zero_HVFs, ]
  )
  Lin_tradeSeq$pseudo <- lineage_subset[[pseudotime_col, drop = TRUE]]
  Lin_tradeSeq$cellWeights <- as.matrix(rep(1, length(Lin_tradeSeq$pseudo)))
  Lin_tradeSeq$genotype <- as.factor(lineage_subset$genotype)

  Lin_tradeSeq$fitGam <- fitGAM(
    counts = Lin_tradeSeq$counts,
    pseudotime = Lin_tradeSeq$pseudo,
    cellWeights = Lin_tradeSeq$cellWeights,
    conditions = as.factor(Lin_tradeSeq$genotype),
    parallel = TRUE
  )
  
  Lin_tradeSeq_diff_exp <- conditionTest(Lin_tradeSeq$fitGam) %>%
    dplyr::filter(pvalue <= 0.05)

  saveRDS(
    Lin_tradeSeq,
    file.path(outdir, paste0(lineage_name, "_tradeSeq_diff_exp.rds"))
  )
  
  write.csv(
    Lin_tradeSeq_diff_exp,
    file.path(outdir, paste0(lineage_name, "_tradeSeq_diff_exp_sig.csv")),
    row.names = FALSE
  )
}

run_tradeseq(FL, "ranked_AC_slingshot_pseudotime", "FL_AC", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/")
run_tradeseq(FL, "ranked_ZC_slingshot_pseudotime", "FL_ZC", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/")
run_tradeseq(FL, "ranked_SC_slingshot_pseudotime", "FL_SC", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/")

run_tradeseq(HL, "ranked_AC_slingshot_pseudotime", "HL_AC", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/")
run_tradeseq(HL, "ranked_ZC_slingshot_pseudotime", "HL_ZC", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/")
run_tradeseq(HL, "ranked_SC_slingshot_pseudotime", "HL_SC", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/")

run_tradeseq(PA_integrated, "ranked_PA1_C_slingshot_pseudotime", "PA1_C", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/")
run_tradeseq(PA_integrated, "ranked_PA1_O_slingshot_pseudotime", "PA1_O", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/")
run_tradeseq(PA_integrated, "ranked_PA2_C_slingshot_pseudotime", "PA2_C", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/")
