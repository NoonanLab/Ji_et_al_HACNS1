library(Seurat)
library(SeuratWrappers)   
library(monocle3)
library(Matrix)
library(dplyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
set.seed(123) 

run_monocle3_modules <- function(seu, pt_col, output_csv) {
  
  cds <- as.cell_data_set(seu)
  rowData(cds)$gene_short_name <- rownames(cds)
  
  cds@colData <- cbind(
    colData(cds),
    seu@meta.data[, setdiff(colnames(seu@meta.data), colnames(colData(cds))), drop = FALSE]
  )
  
  colData(cds)$pt <- seu@meta.data[[pt_col]]
  
  cds <- preprocess_cds(cds, num_dim = 50, method = "PCA")
  cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = TRUE)
  
  hvf <- VariableFeatures(seu)
  is_ribo <- grepl("^(Rps|Rpl|Mrps|Mrpl)", hvf, perl = TRUE)
  is_hsp  <- grepl("^(Hsp|Hspa|Hspb|Hspc|Hspd|Hspe|Hsph|Dnaja|Dnajb|Dnajc)", hvf, perl = TRUE)
  is_mt   <- grepl("^(MT-|mt-)", hvf, perl = TRUE)
  
  hvf <- hvf[!(is_ribo | is_hsp | is_mt)]
  
  SCT <- GetAssayData(seu, assay = "SCT", slot = "data")
  genes_in_SCT <- intersect(hvf, rownames(SCT))
  
  m_sub <- SCT[genes_in_SCT, , drop = FALSE]
  keep_min_gt1 <- if (inherits(m_sub, "Matrix")) {
    Matrix::rowSums(m_sub > 1, na.rm = TRUE) >= 30
  } else {
    rowSums(m_sub > 1, na.rm = TRUE) >= 30
  }
  
  keep_genes <- genes_in_magic[keep_min_gt1]
  cds <- cds[keep_genes, ]
  
  cds <- estimate_size_factors(cds)
  
  fit <- fit_models(
    cds,
    model_formula_str = "~ splines::ns(pt, df = 3)",
    clean_model = TRUE,
    cores = 1
  )
  
  coefs <- coefficient_table(fit)
  
  pt_sig <- coefs %>%
    filter(term != "(Intercept)") %>%
    group_by(gene_id) %>%
    summarize(
      gene_short_name = dplyr::first(gene_short_name),
      min_q_value = suppressWarnings(min(q_value, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(min_q_value)
  
  pt_genes <- pt_sig %>%
    filter(is.finite(min_q_value), min_q_value < 0.05)
  
  cat("Trajectory-variable genes (q < 0.05):", nrow(pt_genes), "\n")
  
  modules <- find_gene_modules(
    cds[pt_genes$gene_id, ],
    resolution = 1e-2,
    random_seed = 123,
    reduction_method = "UMAP",
    cores = 1
  )
  
  modules <- modules %>%
    left_join(
      as_tibble(rowData(cds), rownames = "id")[, c("id", "gene_short_name")],
      by = "id"
    ) %>%
    rename(gene_id = id) %>%
    arrange(module, gene_short_name)
  
  write.csv(modules, output_csv, row.names = FALSE)
  
  modules_tbl <- modules %>%
    transmute(gene_id, gene = gene_short_name, module)
  
  agg <- aggregate_gene_expression(
    cds,
    gene_group_df = modules_tbl %>% dplyr::select(gene_id, module),
    norm_method = "size_only"
  )
  
  module_expr <- if (is.list(agg) && "norm_exp" %in% names(agg)) agg$norm_exp else as.matrix(agg)
  
  ord <- order(colData(cds)$pt)
  pt_vec <- colData(cds)$pt[ord]
  module_expr <- module_expr[, ord, drop = FALSE]
  
  bins <- cut(pt_vec, breaks = 100, labels = FALSE, include.lowest = TRUE)
  
  module_bin_means <- vapply(
    split(seq_along(pt_vec), bins),
    function(idx) rowMeans(module_expr[, idx, drop = FALSE]),
    FUN.VALUE = numeric(nrow(module_expr))
  )
  colnames(module_bin_means) <- paste0("bin", seq_len(ncol(module_bin_means)))
  
  M <- t(scale(t(module_bin_means)))
  M[is.na(M)] <- 0
  
  col_fun <- circlize::colorRamp2(
    seq(-2, 2, length.out = 11),
    RColorBrewer::brewer.pal(11, "RdBu")
  )
  
  Heatmap(
    M,
    name = "z-score",
    col = col_fun,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_title = "Modules",
    column_title = "Pseudotime (binned)"
  )
  
  return(list(
    cds = cds,
    fit = fit,
    coefs = coefs,
    pt_sig = pt_sig,
    pt_genes = pt_genes,
    modules = modules,
    module_expr = module_expr,
    heatmap_matrix = M
  ))
}

FL_Cartilage <- readRDS("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/FL_Cartilage.rds")
DefaultAssay(FL_Cartilage) <- "SCT"
AC_lineage <- FL_Cartilage[ ,!is.na(FL_Cartilage$AC_lineage)]
ZC_lineage <- FL_Cartilage[ ,!is.na(FL_Cartilage$ZC_lineage)]
SC_lineage <- FL_Cartilage[ ,!is.na(FL_Cartilage$SC_lineage)]

run_monocle3_modules(AC_lineage, "ranked_AC_slingshot_pseudotime", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/AC_monocle3_modules.csv")
run_monocle3_modules(ZC_lineage, "ranked_ZC_slingshot_pseudotime", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/ZC_monocle3_modules.csv")
run_monocle3_modules(SC_lineage, "ranked_SC_slingshot_pseudotime", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/FL/SC_monocle3_modules.csv")

HL_Cartilage <- readRDS("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/HL_Cartilage.rds")
DefaultAssay(HL_Cartilage) <- "SCT"
AC_lineage <- HL_Cartilage[ ,!is.na(HL_Cartilage$AC_lineage)]
ZC_lineage <- HL_Cartilage[ ,!is.na(HL_Cartilage$ZC_lineage)]
SC_lineage <- HL_Cartilage[ ,!is.na(HL_Cartilage$SC_lineage)]

run_monocle3_modules(AC_lineage, "ranked_AC_slingshot_pseudotime", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/AC_monocle3_modules.csv")
run_monocle3_modules(ZC_lineage, "ranked_ZC_slingshot_pseudotime", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/ZC_monocle3_modules.csv")
run_monocle3_modules(SC_lineage, "ranked_SC_slingshot_pseudotime", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/HL/SC_monocle3_modules.csv")

PA_Osteoblasts <- readRDS("/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA_Osteoblasts.rds")
PA1_C_lineage <- PA_Osteoblasts[ ,!is.na(PA_Osteoblasts$PA1_C_lineage)]
PA1_O_lineage <- PA_Osteoblasts[ ,!is.na(PA_Osteoblasts$PA1_O_lineage)]
PA2_C_lineage <- PA_Osteoblasts[ ,!is.na(PA_Osteoblasts$PA2_C_lineage)]
PA2_C_lineage <- PA2_C_lineage[ ,PA2_C_lineage$stage %in% c("E9.5", "E10.5", "E11.5")]

run_monocle3_modules(PA1_C_lineage, "ranked_PA1_C_slingshot_pseudotime", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA1_C_monocle3_modules.csv")
run_monocle3_modules(PA1_O_lineage, "ranked_PA1_O_slingshot_pseudotime", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA1_O_monocle3_modules.csv")
run_monocle3_modules(PA2_C_lineage, "ranked_PA2_C_slingshot_pseudotime", "/vast/palmer/scratch/noonan/yj345/new/CCR/integrated/Tissue/PA/PA2_C_monocle3_modules.csv")
