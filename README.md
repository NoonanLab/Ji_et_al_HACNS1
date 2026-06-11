Ji et al. - A Human Accelerated Region Drives Opposing Heterochronic Changes in Craniofacial and Limb Development

ATAC-seq
Step 1
ModifyGenome.py and Bowtie2Index.sh – Generate genotype-specific Bowtie2 indices by replacing the corresponding mm10 locus with either the human HACNS1 sequence or the chimpanzee ortholog sequence.
Step 2
ATAC.sh – Perform quality control of ATAC-seq reads using FastQC and trim adapter sequences using cutadapt. Align ATAC-seq reads to the appropriate genotype-specific genome using Bowtie2. Sort aligned reads using SAMtools, remove duplicate reads using Picard, exclude mitochondrial reads prior to peak calling, and call ATAC-seq narrow peaks using MACS2.
Step 3
DiffBind.R – Perform differential chromatin accessibility analysis using DiffBind.
snRNA-seq
Step 4
run_barcode.sh – Generate barcode whitelists for each snRNA-seq dataset using PIPseeker barcode.
Step 5
run_velocyto.sh – Align snRNA-seq reads to the mm10 reference genome using STAR.
Step 6
CreateSeuratObject.R, CellCycleRegression.R, and DoubletFinder.R – Perform initial snRNA-seq quality control in Seurat, including filtering cells by detected genes and mitochondrial content, cell-cycle scoring, doublet removal using DoubletFinder, and SCTransform normalization.
Step 7
Integrate.R – Integrate samples across developmental stages using Harmony.
Step 8
3d_phate.py – Generate three-dimensional PHATE embeddings to visualize global developmental structure in integrated mesenchymal cells.
Step 9
URD.R – Reconstruct developmental trajectories within integrated mesenchymal populations using URD diffusion maps and biased random walks from terminal cell states to defined root populations.
Step 10
Phate_after_URD.py and Slingshot_after_URD – Subset cells with high random-walk visitation frequencies, generate two-dimensional PHATE embeddings, and infer lineage trajectories and pseudotime using pyslingshot.
Step 11
Monocle3.R – Identify genes with trajectory-dependent expression patterns using the Monocle3 graph_test function and cluster significant genes into modules using Louvain clustering.
Step 12
cellalign.R – Align developmental trajectories between genotypes using CellAlign.
Step 13
tradeseq.R – Identify genotype-dependent transcriptional shifts along developmental trajectories using tradeSeq.
Step 14
ClusterGenes.R – Cluster significantly shifted genes.
Step 15
GO.R – Perform Gene Ontology enrichment analysis.
Gene regulatory network construction
Step 16
run_pipeline.sh – Re-align E10.5 ATAC-seq data to the mm10 genome before footprint analysis.
Step 17
ChIPseeker.R – Annotate promoter-proximal ATAC-seq peaks using ChIPseeker.
Step 18
DownloadCapHiC.sh and DownloadChIPSeqData.sh – Download published Capture Hi-C and E10.5 H3K27ac ChIP-seq datasets.
Step 19
EnhancerAnnotation.sh and preCreateNetWork.sh – Annotate active enhancers using ATAC-seq, H3K27ac ChIP-seq, and Capture Hi-C interactions.
Step 20
TOBIAS.sh – Perform footprint analysis using TOBIAS. First, correct ATAC-seq signal for Tn5 insertion bias using TOBIAS-ATACorrect and compute base-pair footprint scores using TOBIAS-FootprintScores. Second, identify transcription factor binding motifs and genotype-specific differences in putative transcription factor binding using TOBIAS-BINDetect. Third, construct transcription factor binding networks using annotated ATAC-seq peaks and transcription factor binding sites with TOBIAS CreateNetwork.
Step 21
Magic.py – Denoise snRNA-seq expression matrices using MAGIC prior to coexpression-based gene regulatory network inference.
Step 22
DownSampling.py, runGRN.sh, and CombineGRNResults.py – Infer transcription factor-target coexpression modules using GRNBoost2 implemented in pySCENIC.
Step 23
Combine_grn_TOBIAS.R – Filter direct binding networks using high-confidence coexpression modules to retain transcription factor-target interactions.
Step 24
PermutationTest.R – Test whether genotype-shifted genes are enriched at specific downstream distances from Gbx2 in tissue-specific gene regulatory networks using two-tailed permutation testing.

