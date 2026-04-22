library(Seurat)
library(tradeSeq)
library(cellAlign)
library(dplyr)
library(Matrix)
library(pheatmap)

run_cluster_by_module <- function(
  rds_file,
  csv_file,
  outdir,
  k_req = 2
) {
  x <- readRDS(rds_file)
  interScaledGlobalHumRTN <- x$Human
  interScaledGlobalChRTN  <- x$Chimp
  
  df <- read.csv(csv_file, row.names = 1, check.names = FALSE)
  
  hum  <- interScaledGlobalHumRTN$scaledData
  chim <- interScaledGlobalChRTN$scaledData
  
  common_genes <- intersect(rownames(hum), rownames(chim))
  hum  <- hum[common_genes, , drop = FALSE]
  chim <- chim[common_genes, , drop = FALSE]
  
  df$kmeans_cluster <- NA_integer_
  mods <- sort(unique(df$module))
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (m in mods) {
    genes_m <- intersect(rownames(df)[df$module == m], common_genes)
    
    # drop genes with any NA across the interpolated curves
    if (length(genes_m) > 0) {
      keep <- rowSums(is.na(hum[genes_m, , drop = FALSE])) == 0 &
        rowSums(is.na(chim[genes_m, , drop = FALSE])) == 0
      genes_m <- genes_m[keep]
    }
    
    # need at least k_req genes
    if (length(genes_m) < k_req) {
      message(sprintf("Module %s: %d gene(s) < k=%d — skipping.", m, length(genes_m), k_req))
      next
    }
    
    res <- tryCatch(
      pseudotimeClust(
        x = hum[genes_m, , drop = FALSE],
        y = chim[genes_m, , drop = FALSE],
        k = k_req
      ),
      error = function(e) e
    )
    
    if (inherits(res, "error")) {
      message(sprintf("Module %s: clustering failed (%s) — skipping.", m, conditionMessage(res)))
      next
    }
    
    cl <- NULL
    if (!is.null(res$cluster)) cl <- res$cluster
    else if (!is.null(res$clustering)) cl <- res$clustering
    else if (!is.null(res$kmeansRes$cluster)) cl <- res$kmeansRes$cluster
    
    if (is.null(cl)) {
      message(sprintf("Module %s: no cluster labels found — skipping.", m))
      next
    }
    
    if (is.null(names(cl))) names(cl) <- genes_m
    cl <- cl[genes_m]
    df[genes_m, "kmeans_cluster"] <- as.integer(cl)
    
    uniq_k <- sort(unique(as.integer(cl)))
    for (kk in uniq_k) {
      obj_to_save <- list(
        module         = m,
        kmeans_cluster = kk,
        genes          = names(cl)[cl == kk],
        center         = if (!is.null(res$centers) && is.matrix(res$centers)) res$centers[kk, , drop = TRUE] else NULL,
        kmeans_result  = res
      )
      
      fpath <- file.path(outdir, sprintf("module_%s_kmeans_%d.rds", m, kk))
      saveRDS(obj_to_save, fpath)
    }
  }
  
  write.csv(df, csv_file)
}

run_cluster_by_module(
  rds_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_AC_interScaledGlobal_RTNs.rds",
  csv_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_AC_tradeSeq_diff_exp_l2fc_58_sig.csv",
  outdir   = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_AC/kmeans_cluster_by_module"
)
run_cluster_by_module(
  rds_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_ZC_interScaledGlobal_RTNs.rds",
  csv_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_ZC_tradeSeq_diff_exp_l2fc_58_sig.csv",
  outdir   = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_ZC/kmeans_cluster_by_module"
)
run_cluster_by_module(
  rds_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_SC_interScaledGlobal_RTNs.rds",
  csv_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_SC_tradeSeq_diff_exp_l2fc_58_sig.csv",
  outdir   = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_SC/kmeans_cluster_by_module"
)

run_cluster_by_module(
  rds_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_AC_interScaledGlobal_RTNs.rds",
  csv_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_AC_tradeSeq_diff_exp_l2fc_58_sig.csv",
  outdir   = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_AC/kmeans_cluster_by_module"
)
run_cluster_by_module(
  rds_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_ZC_interScaledGlobal_RTNs.rds",
  csv_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_ZC_tradeSeq_diff_exp_l2fc_58_sig.csv",
  outdir   = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_ZC/kmeans_cluster_by_module"
)
run_cluster_by_module(
  rds_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_SC_interScaledGlobal_RTNs.rds",
  csv_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_SC_tradeSeq_diff_exp_l2fc_58_sig.csv",
  outdir   = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_SC/kmeans_cluster_by_module"
)

run_cluster_by_module(
  rds_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA1_C_interScaledGlobal_RTNs.rds",
  csv_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA1_C_tradeSeq_diff_exp_l2fc_58_sig.csv",
  outdir   = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA1_C/kmeans_cluster_by_module"
)
run_cluster_by_module(
  rds_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA1_O_interScaledGlobal_RTNs.rds",
  csv_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA1_O_tradeSeq_diff_exp_l2fc_58_sig.csv",
  outdir   = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA1_O/kmeans_cluster_by_module"
)
run_cluster_by_module(
  rds_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA2_C_interScaledGlobal_RTNs.rds",
  csv_file = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA2_C_tradeSeq_diff_exp_l2fc_58_sig.csv",
  outdir   = "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA2_C/kmeans_cluster_by_module"
)
