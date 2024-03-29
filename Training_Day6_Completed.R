### Training Day 6 ###
## By Cole Nawrocki ##

## Topic: Practice Assignment ## 

# Data from this paper: 
# https://www.nature.com/articles/s41467-022-31949-2#data-availability

# Data to Use: 
# Look in ~/Bioinformatics/Learning-Resources/Example-Data/day6_data

# Packages to Use: 
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(tidyverse)
library(Matrix)
library(DESeq2)
library(data.table)
library(scDblFinder)
library(SeuratDisk)
library(ggrepel)

# Start with reading the data. HINT: use fread() function from the data.table
# package. 
data_dir <- "/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Grainger-Lab/Bioinformatics/Learning-Resources/Example-Data/day6_data"
files <- list.files(data_dir)
all_cts <- list()
for (f in files[grep(".txt", files)]) {
  cts <- fread(paste(data_dir, f, sep = "/"))
  cts <- column_to_rownames(cts, "GENE")
  all_cts[[gsub("_dge.txt","",f)]] <- as(as.matrix(cts), "CsparseMatrix")
}

meta <- read.csv(paste(data_dir, "St48_cell_info.csv", sep = "/"), row.names = 1)

# Make sure all of the cell IDs match in the metadata and in the counts table.
all_meta <- list()
for (nm in names(all_cts)) {
  colnames(all_cts[[nm]]) <- paste(nm, colnames(all_cts[[nm]]), sep = ".")
  all_meta[[nm]] <- meta[grep(nm, rownames(meta), fixed = T),]
}

for (nm in names(all_cts)) {
  all_cts[[nm]] <- all_cts[[nm]][,rownames(all_meta[[nm]])]
}

for (nm in names(all_cts)) {
  print(all(rownames(all_meta[[nm]]) == colnames(all_cts[[nm]])))
}

# Next, do QC for each sample individually. 
objs <- list()
for (nm in names(all_cts)) {
  objs[[nm]] <- CreateSeuratObject(counts = all_cts[[nm]], 
                                   meta.data = all_meta[[nm]], 
                                   project = "xenopus_training")
}

pp_qc <- function(obj_name, mad_cutoff = 3) {
  
  seurat_obj <- objs[[obj_name]]
  
  seurat_obj$pct_mt <- PercentageFeatureSet(seurat_obj, pattern = "^mt")
  top20 <- names(sort(rowSums(seurat_obj@assays$RNA@counts), decreasing = T))[1:20]
  seurat_obj$pct_counts_inTop20_genes <- PercentageFeatureSet(seurat_obj, features = top20)
  seurat_obj$log1p_nCount_RNA <- log1p(seurat_obj$nCount_RNA)
  seurat_obj$log1p_nFeature_RNA <- log1p(seurat_obj$nFeature_RNA)
  
  outlier_df <- data.frame(row.names = rownames(seurat_obj@meta.data))
  for (metric in c("pct_counts_inTop20_genes", "log1p_nCount_RNA", "log1p_nFeature_RNA", "pct_mt")) {
    M <- seurat_obj@meta.data[,c(metric)]
    names(M) <- rownames(seurat_obj@meta.data)
    mean_abs_dev <- mad(M)
    outlier <- (M < median(M)-mad_cutoff*mean_abs_dev) | (M > median(M)+mad_cutoff*mean_abs_dev)
    outlier_df[[metric]] <- outlier
  }
  
  doublets <- scDblFinder(seurat_obj@assays$RNA@counts)
  seurat_obj <- AddMetaData(seurat_obj, metadata = doublets@colData$scDblFinder.class, col.name = "multiplet_class")
  
  p1 <- VlnPlot(seurat_obj, features = c("pct_counts_inTop20_genes", "nCount_RNA", "nFeature_RNA", "pct_mt"), ncol = 4)
  seurat_obj <- AddMetaData(seurat_obj, rowSums(outlier_df) > 0, col.name = "outlier_status")
  seurat_obj <- subset(seurat_obj, outlier_status == FALSE)
  p2 <- VlnPlot(seurat_obj, features = c("pct_counts_inTop20_genes", "nCount_RNA", "nFeature_RNA", "pct_mt"), ncol = 4)
  
  return(list(seurat_obj, p1, p2))
}

pp_qc_result <- lapply(names(objs), pp_qc)
names(pp_qc_result) <- names(objs)

for (samp in names(pp_qc_result)) {
  cat(samp, "\n", sep = "")
  print(pp_qc_result[[samp]][[2]])
  ggsave(paste(samp, "before_qc.png", sep = "_"), width = 16, height = 8, units = "in", bg="white")
  print(pp_qc_result[[samp]][[3]])
  ggsave(paste(samp, "after_qc.png", sep = "_"), width = 16, height = 8, units = "in", bg="white")
}

# Combine the samples. 
st48 <- merge(x = pp_qc_result[[1]][[1]], 
              y = c(pp_qc_result[[2]][[1]], pp_qc_result[[3]][[1]], pp_qc_result[[4]][[1]]),
              project = "xenopus_training")

# Filter out genes that do not appear in at least 10 cells.
useable_genes <- which((st48@assays$RNA@counts != 0) |> rowSums() > 10)
st48 <- subset(st48, features = useable_genes)

# Seurat standard workflow: Normalizing and scaling, clustering, etc. 
st48 <- NormalizeData(st48, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
st48 <- FindVariableFeatures(st48, selection.method = "vst", nfeatures = 3000, verbose = F)
FeatureScatter(st48, feature2 = "nCount_RNA", feature1 = "pct_mt", group.by = "orig.ident")
FeatureScatter(st48, feature2 = "nCount_RNA", feature1 = "pct_counts_inTop20_genes", group.by = "orig.ident")
options(future.globals.maxSize = 4e9)
plan(strategy = "multisession", workers=8)
st48 <- ScaleData(st48, verbose = F)
st48 <- RunPCA(st48, features = VariableFeatures(st48), npcs = 60, verbose = F)
plan(strategy = "sequential")
ElbowPlot(st48, ndims = 60)
st48 <- RunUMAP(st48, dims = 1:30, verbose = F)

DimPlot(st48, group.by = "multiplet_class")

st48$sample_id <- substr(st48$cellID, 1, 12)
DimPlot(st48, group.by = "sample_id")
DimPlot(st48, group.by = "celltype")

st48 <- FindNeighbors(st48, dims = 1:30, verbose = F)
st48 <- FindClusters(st48, algorithm = 1, resolution = seq(0.1, 3, 0.1), verbose = F)

for (i in seq(0.1, 3, 0.1)) {
  print(DimPlot(st48, group.by = paste("RNA_snn_res",i, sep = ".")) + NoLegend())
}

st48$clusters_44 <- st48$RNA_snn_res.1.5
st48 <- SetIdent(st48, value = "clusters_44")
DimPlot(st48, label = T) + NoLegend()

# Save the data as .Rds and as .h5Seurat
SaveH5Seurat(st48, "xenopus_training.h5Seurat")
SaveSeuratRds(st48, "xenopus_training.Rds")

# Reading the data in (I am picking back up here at a later date).
st48 <- LoadSeuratRds("xenopus_training.Rds")

# Differential expression analysis to find cluster markers. 
cluster_markers <- FindAllMarkers(st48, logfc.threshold = 1, test.use = "wilcox", only.pos = T)
top_10s <- cluster_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC)

# Differential expression analysis with pseudobulking to compare clusters 10 and 13. 
# Subseting for the cells we want.
cl10_cl13 <- subset(st48, clusters_44 %in% c(10,13))

# Checking it worked.
DimPlot(cl10_cl13)

# Aggregating the counts.
psb_cts <- AggregateExpression(cl10_cl13, group.by = c("clusters_44", "sample_id"))[["RNA"]]

# What the counts data looks like now.
psb_cts[1:4,1:4]

# Forming the new metadata. 
psb_meta <- data.frame(row.names = colnames(psb_cts))
psb_meta$cluster <- str_split(rownames(psb_meta), pattern = "_", simplify = T)[,1]
psb_meta$cluster <- ifelse(psb_meta$cluster == "g10", yes = "cluster10", no = "cluster13")

# What the metadata looks like now. 
psb_meta

# Running DESeq2. 
dds <- DESeqDataSetFromMatrix(countData = psb_cts, 
                              colData = psb_meta, 
                              design = ~cluster)
dds <- DESeq(dds)

# Obtaining the results table. Note: the contrast argument is very useful. Here, 
# I have specified that positive LFC will correspond to cluster 10 and negative 
# LFC will correspond to cluster 13.
res <- results(dds, contrast = c("cluster", "cluster10", "cluster13")) |> as.data.frame() |> na.omit()

# Summarizing with a volcano plot. 
res$delabel<-'Neither'
res$delabel[res$padj<0.05 & res$log2FoldChange > 1]<-"cluster10"
res$delabel[res$padj<0.05 & res$log2FoldChange < -1]<-"cluster13"
res <- res[order(res$delabel),]
res$lbl<-NA
res$lbl[res$delabel != 'Neither'] <- rownames(res)[res$delabel != 'Neither']
v.plot<-ggplot(data=res,aes(x=log2FoldChange, y=-log10(padj),col=delabel, label=lbl)) +
  geom_point() + theme_minimal() + geom_text_repel(max.overlaps = 50) +
  labs(title="Cluster 13 vs Cluster 10") + 
  scale_color_manual(values=c("cluster10"="red", "cluster13"="blue","Neither"="gray")) +
  geom_vline(xintercept=c(-1,1),col="black") + 
  geom_hline(yintercept=-log10(0.05),col="black")
print(v.plot)

# Looking at a gene that was identified as up-regulated in cluster 13 compared 
# to cluster 10. 
FeaturePlot(cl10_cl13, features = "vcam1.S", label = T, order = T)
plotCounts(dds, gene = "vcam1.S", intgroup = "cluster")
