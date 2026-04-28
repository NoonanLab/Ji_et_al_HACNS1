library(DiffBind)

samples <- read.csv("path/to/samples.csv")
dba_object <- dba(sampleSheet = samples)

dba.plotHeatmap(dba_object)
dba.plotPCA(dba_object, label = DBA_TISSUE)

dba_object <- dba.count(dba_object, summits = 250)
dba_object <- dba.contrast(dba_object, categories = DBA_TISSUE)

dba_object <- dba.analyze(dba_object)

diff_results <- dba.report(dba_object)
write.csv(as.data.frame(diff_results), file = "diffbind_results.csv")
