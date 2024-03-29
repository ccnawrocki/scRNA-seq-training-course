### Training Day 2 ###
## By Cole Nawrocki ##

## Topic: Seurat Workflow Part 2 -- Adding QC and finishing the process ##

# Load packages. To install scDblFinder, see Training_Day3.R, which is on 
# packages and package managers. 
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(scDblFinder)

# Set the objects to be the older, simpler version.
options(Seurat.object.assay.version = "v3")

# Read in data
cts <- Read10X("/Users/anwita.molaka/Library/CloudStorage/OneDrive-UniversityofVirginia/Bioinformatics/Learning-Resources/Example-Data/filtered_feature_bc_matrix")
cts <- cts[[1]]

# Last time, we didn't do QC. Normally, you want to. 
# You want to filter out outlier cells and genes. You also want to detect possible doublets (two cells mistaken as one). 
# For a long time, people would just make the Seurat object, look at violin plots and set arbitrary cutoffs. 
# Now people use the MAD statistic to set the cutoffs. 
obj <- CreateSeuratObject(counts = cts, project = "Training2")

# Looking at some QC features. People used to just eyeball it and cut it at, say, 15000 for nCount_RNA and 4000 for nFeature_RNA.
VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA"))

# First, lets calculate the necessary metrics for our new methods. I got these from Single Cell Best Practices. 
obj$pct_mt <- PercentageFeatureSet(obj, pattern = "^MT")
top20 <- names(sort(rowSums(obj@assays$RNA@counts), decreasing = T))[1:20]
obj$pct_counts_inTop20_genes <- PercentageFeatureSet(obj, features = top20)
obj$log1p_nCount_RNA <- log1p(obj$nCount_RNA)
obj$log1p_nFeature_RNA <- log1p(obj$nFeature_RNA)

# Identifying possible outliers. We are using MAD > 5. 
# A note on for loops. If you aren't looping through a ton of things, they are fine. Don't listen to the STAT department. 
outlier_df <- data.frame(row.names = rownames(obj@meta.data))
for (metric in c("pct_counts_inTop20_genes", "log1p_nCount_RNA", "log1p_nFeature_RNA", "pct_mt")) {
  M <- obj@meta.data[,c(metric)]
  names(M) <- rownames(obj@meta.data)
  mean_abs_dev <- mad(M)
  outlier <- (M < median(M)-5*mean_abs_dev) | (M > median(M)+5*mean_abs_dev)
  print(sum(outlier))
  outlier_df[[metric]] <- outlier
}

# Filtering out the outliers. 
obj <- AddMetaData(obj, rowSums(outlier_df) > 0, col.name = "outlier_status")
obj <- subset(obj, outlier_status == FALSE)
VlnPlot(obj, features = c("pct_counts_inTop20_genes", "nCount_RNA", "nFeature_RNA", "pct_mt"), ncol = 4)

# Doing doublet detection. We will decide if they need to be removed later. 
doublets <- scDblFinder(obj@assays$RNA@counts)
obj <- AddMetaData(obj, metadata = doublets@colData$scDblFinder.class, col.name = "multiplet_class")

# Same steps as last week.
obj <- NormalizeData(obj, normalization.method = "LogNormalize")
obj <- FindVariableFeatures(obj, method = "vst", nfeatures = 2000)
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "pct_mt")
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "pct_counts_inTop20_genes")

options(future.globals.maxSize = 2e9)
plan(strategy = "multisession", workers=4)

obj <- ScaleData(obj, vars.to.regress = "nFeature_RNA")
obj <- RunPCA(obj, npcs = 50)

ElbowPlot(obj, ndims = 50)

obj <- RunUMAP(obj, dims = 1:20)

plan(strategy = "sequential")

# Should we remove doublets? They all sort of cluster together, so YES. 
DimPlot(obj, group.by = "multiplet_class")

# Remove them and redo the previous steps. 
obj <- subset(obj, multiplet_class == "singlet")
obj <- NormalizeData(obj, normalization.method = "LogNormalize")
obj <- FindVariableFeatures(obj, method = "vst", nfeatures = 2000)
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "pct_mt")
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "pct_counts_inTop20_genes")
options(future.globals.maxSize = 2e9)
plan(strategy = "multisession", workers=4)
obj <- ScaleData(obj, vars.to.regress = "nFeature_RNA")
obj <- RunPCA(obj, npcs = 50)
ElbowPlot(obj, ndims = 50)
obj <- RunUMAP(obj, dims = 1:20)
plan(strategy = "sequential")

# Now the data looks like this. 
DimPlot(obj, group.by = "multiplet_class")

# New steps. 
# Finding neighbors. 
obj <- FindNeighbors(obj, dims = 1:20)

# Clustering. Can do this iteratively. We are using the "Louvain" algorithm. The "Leiden" algorithm is good too. 
obj <- FindClusters(obj, algorithm = 1, resolution = seq(0.5,3,0.5))
for (i in seq(0.5,3,0.5)) {
  print(DimPlot(obj, group.by = paste("RNA_snn_res", i, sep = "."), label = T))
}

# Resolution of 1 looks good. This is sort of an arbitrary value though. 
# What I find best is to look for a marker gene using FeaturePlot. 
# Pick the resolution that confines the marker gene to a cluster.
obj$seurat_clusters <- obj$RNA_snn_res.1
ob <- SetIdent(obj, value = "seurat_clusters")
DimPlot(obj)

# Finding markers. We will discuss the arguments at a later date. 
de <- FindAllMarkers(obj, logfc.threshold = 1, test.use = "wilcox", only.pos = T)
FeaturePlot(obj, "CD14", label = T)

# Saving the data. 
# As .h5Seurat (Need SeuratDisk package, but more reliable, apparently)
SaveH5Seurat(obj, "/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Grainger-Lab/Bioinformatics/Learning-Resources/Example-Data/Training2.h5Seurat")

# As .Rds (Don't need SeuratDisk package)
SaveSeuratRds(obj, "/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Grainger-Lab/Bioinformatics/Learning-Resources/Example-Data/Training2.Rds")

# As .h5ad (AnnData format for opening in Scanpy)
SeuratDisk::Convert("/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Grainger-Lab/Bioinformatics/Learning-Resources/Example-Data/Training2.h5Seurat", 
                    "/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Grainger-Lab/Bioinformatics/Learning-Resources/Example-Data/Training2.h5ad")

# Loading the data. 
# From .h5Seurat
obj_test1 <- LoadH5Seurat("/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Grainger-Lab/Bioinformatics/Learning-Resources/Example-Data/Training2.h5Seurat")

# From .Rds
obj_test2 <- LoadSeuratRds("/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Grainger-Lab/Bioinformatics/Learning-Resources/Example-Data/Training2.Rds")

