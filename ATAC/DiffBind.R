library(DiffBind)

all_samples_file <- "~/palmer_scratch/ATAC/atac_samples.csv"
outdir <- "~/palmer_scratch/ATAC/DiffBind_all_samples"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

all_samples <- read.csv(all_samples_file, stringsAsFactors = FALSE)

run_diffbind_one_group <- function(tissue, timepoint) {
  prefix <- paste0(tissue, "_", timepoint)

  samples <- subset(
    all_samples,
    Tissue == tissue & Timepoint == timepoint
  )

  sample_outfile <- file.path(outdir, paste0("sample_", prefix, ".csv"))
  write.csv(samples, sample_outfile, row.names = FALSE, quote = FALSE)

  dbobj_file <- file.path(outdir, paste0("dbObj_DiffBind_", prefix, ".rds"))
  report_file <- file.path(outdir, paste0("dbObj_DB_DiffBind_", prefix, ".rds"))
  report_csv <- file.path(outdir, paste0("DiffBind_report_", prefix, ".csv"))
  plot_pdf <- file.path(outdir, paste0("DiffBind_plot_", prefix, ".pdf"))

  dbObj <- dba(sampleSheet = samples, scoreCol = 5)

  print(dbObj)

  pdf(plot_pdf)
  plot(dbObj)
  dev.off()

  dbObj <- dba.count(dbObj, bParallel = TRUE)
  saveRDS(dbObj, file = dbobj_file)

  dbObj <- dba.normalize(dbObj)
  saveRDS(dbObj, file = dbobj_file)

  dbObj <- dba.contrast(
    dbObj,
    contrast = c("Condition", "Human", "Chimp"),
    minMembers = 2
  )

  dbObj <- dba.analyze(dbObj, bParallel = TRUE)

  saveRDS(dbObj, file = dbobj_file)

  dbObj_DB <- dba.report(dbObj)

  saveRDS(dbObj_DB, file = report_file)

  dbObj_DB_df <- as.data.frame(dbObj_DB)
  write.csv(dbObj_DB_df, report_csv, row.names = FALSE)

  return(dbObj)
}

groups <- unique(all_samples[, c("Tissue", "Timepoint")])
groups <- groups[order(groups$Tissue, groups$Timepoint), ]

print(groups)

results <- list()

for (i in seq_len(nrow(groups))) {

  tissue <- groups$Tissue[i]
  timepoint <- groups$Timepoint[i]
  prefix <- paste0(tissue, "_", timepoint)

  results[[prefix]] <- run_diffbind_one_group(tissue, timepoint)
}
