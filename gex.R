mdir <- "~/work/partis/datascripts/meta/goo-dengue-10x"
plotdir <- "/fh/fast/matsen_e/dralph/partis/tmp/gex" #"~/bDropbox/tmp-plots"
feature.matrix.fname <- "/fh/fast/matsen_e/data/goo-dengue-10x/data/gex/filtered_feature_bc_matrix"

# install.packages("BiocManager")
## BiocManager::install(c("scRNAseq", "scater", "scran", "uwot"), dependencies=TRUE)  # "TENxPBMCData"
## BiocManager::install("DropletUtils")
library(DropletUtils, warn.conflicts=F, quietly=T)
library(scater, warn.conflicts=F, quietly=T)
library(scran, warn.conflicts=F, quietly=T)
library(pheatmap, warn.conflicts=F, quietly=T)
sce <- read10xCounts(feature.matrix.fname)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

# Quality control.
is.mito <- grepl("^MT-", rownames(sce))  # figure out which genes are mitochondrial
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")  # identifies + removes outliers (in several qc metrics)
sce <- sce[, !filtered$discard]

# Normalization.
sce <- logNormCounts(sce)

# Feature selection
# genes from fabio (200 most discriminatory between plasmablast + naive B cell):
## fabio.pb.genes <- read.csv(sprintf("%s/plasmablast_markers.tsv", mdir), sep="\t", header=T)$GeneName # $name
waick.genes <- read.csv(sprintf("%s/waickman-markers.csv", mdir), header=T)$gene
genelist <- waick.genes  # fabio.pb.genes
print(sprintf("  using %d genes: %s", length(genelist), paste(genelist, collapse=" ")))
gene.bools <- rowData(sce)$Symbol %in% genelist  # $ID

# Dimensionality reduction.
set.seed(1)
n.comp <- min(25, as.integer(length(genelist)/2))
print(sprintf("running pca with %d components", n.comp))
sce <- runPCA(sce, ncomponents=n.comp, subset_row=gene.bools)
sce <- runUMAP(sce, dimred='PCA', external_neighbors=TRUE)  # uses pca results from previous step TODO test variety of N neighbors and min_dist values

# Clustering.
g <- buildSNNGraph(sce, use.dimred='PCA')
colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)

# write output files (more written below)
capture.output(reducedDims(sce)[[1]], file=sprintf("%s/pca.txt", plotdir))
capture.output(reducedDims(sce)[[2]], file=sprintf("%s/umap.txt", plotdir))
capture.output(colLabels(sce), file=sprintf("%s/clusters.txt", plotdir))

## pdf(sprintf("%s/clusters.pdf", plotdir))
png(sprintf("%s/clusters.png", plotdir))
plotUMAP(sce, colour_by="label")
dev.off()

# find marker genes
## pval.type="all" looks for genes that distinguish the cluster from *all* clusters, rather than the default of looking in any pairwise comparison
markers <- findMarkers(sce)  # <markers>: list of data frames for each cluster NOTE this uses *all* the genes, and i can't figure out a way to tell it not to
n.genes = 10
print(sprintf("  top %d genes for each cluster (total size %d)", n.genes, length(sce$label)))
for(ich in seq(length(markers))) {  # look at genes that distinguish cluster ich from all other clusters
    print(sprintf("   cluster %2d  size %4d  frac %.2f", ich, sum(sce$label==ich), sum(sce$label==ich) / length(sce$label)))
    interesting <- markers[[ich]]
    capture.output(interesting[1:n.genes,], file=sprintf("%s/markers-cluster-%d.txt", plotdir, ich))
    best.set <- interesting[interesting$Top <= n.genes,]  # look at the top N genes from each pairwise comparison
    logFCs <- getMarkerEffects(best.set)
    png(sprintf("%s/heatmap-%d.png", plotdir, ich))
    pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
    dev.off()
}
