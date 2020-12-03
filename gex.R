# install.packages("BiocManager")
## BiocManager::install(c("scRNAseq", "scater", "scran", "uwot"), dependencies=TRUE)  # "TENxPBMCData"
## BiocManager::install("DropletUtils")
library(DropletUtils, warn.conflicts=F, quietly=T)
library(scater, warn.conflicts=F, quietly=T)
library(scran, warn.conflicts=F, quietly=T)
library(pheatmap, warn.conflicts=F, quietly=T)
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

pdf("~/Dropbox/tmp-plots/clusters.pdf")
plotUMAP(sce, colour_by="label")
dev.off()

# find marker genes
markers <- findMarkers(sce)  # list of data frames for each cluster
n.genes = 3
print(sprintf("  top %d genes for each cluster (total size %d)", n.genes, length(sce$label)))
for(ich in seq(length(markers))) {  # look at genes that distinguish cluster ich from all other clusters
    print(sprintf("   cluster %2d  size %4d  frac %.2f", ich, sum(sce$label==ich), sum(sce$label==ich) / length(sce$label)))
    interesting <- markers[[ich]]
    best.set <- interesting[interesting$Top <= n.genes,]  # look at the top N genes from each pairwise comparison
    logFCs <- getMarkerEffects(best.set)
    png(sprintf("~/Dropbox/tmp-plots/heatmap-%d.png", ich))
    pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
    dev.off()
}
