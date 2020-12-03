# install.packages("BiocManager")
## BiocManager::install(c("scRNAseq", "scater", "scran", "uwot"), dependencies=TRUE)  # "TENxPBMCData"
## BiocManager::install("DropletUtils")
library(DropletUtils, warn.conflicts=F)
library(scater, warn.conflicts=F)
library(scran, warn.conflicts=F)
sce <- read10xCounts("/fh/fast/matsen_e/data/goo-dengue-10x/data/gex/filtered_feature_bc_matrix")
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

# Quality control.
is.mito <- grepl("^MT-", rownames(sce))  # figure out which genes are mitochondrial
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")  # identifies + removes outliers (in several qc metrics)
sce <- sce[, !filtered$discard]

# Normalization.
sce <- logNormCounts(sce)

# Feature selection
## # highly variable genes:
## dec <- modelGeneVar(sce)
## hvg <- getTopHVGs(dec, prop=0.1)
# or genes from fabio (200 most discriminatory between plasmablast + naive B cell):
fgenes <- read.csv("/fh/fast/matsen_e/data/goo-dengue-10x/plasmablast_markers.tsv", sep="\t", header=T)$GeneName # $name
cdf.genes <- rowData(sce)$Symbol %in% fgenes  # $ID

# Dimensionality reduction.
set.seed(1)
sce <- runPCA(sce, ncomponents=25, subset_row=cdf.genes) # =hvg
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)  # uses pca results from previous step TODO test variety of N neighbors and min_dist values

# Clustering.
g <- buildSNNGraph(sce, use.dimred = 'PCA')
colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)

png("~/Dropbox/tmp-plots/clusters.png")
plotUMAP(sce, colour_by="label")
dev.off()

# find marker genes
markers <- findMarkers(sce)  # list of data frames for each cluster
chosen <- "1"  # look at genes that distinguish cluster 1 from all other clusters
interesting <- markers[[chosen]]
best.set <- interesting[interesting$Top <= 6,]  # look at the top 6 genes from each pairwise comparison
logFCs <- getMarkerEffects(best.set)
library(pheatmap, warn.conflicts=F)
png("~/Dropbox/tmp-plots/heatmap.png")
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
dev.off()
