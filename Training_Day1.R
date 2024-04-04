### Training Day 1 ###
## By Cole Nawrocki ##

## Topic: Seurat Workflow Part 1 -- Setting up the object and basic process ##

# Load packages
library(Seurat)
library(Matrix)

# Set the objects to be the older, simpler version.
options(Seurat.object.assay.version = "v3")

# Read in counts. 
# Note: you may have to change the file path accordingly.
cts <- read.csv("/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Grainger-Lab/Bioinformatics/Learning-Resources/Example-Data/Lung12_exprMat_file.csv")
cts

# Get the counts into the correct format.
rownames(cts) <- paste(cts$fov, cts$cell_ID, sep = "_")
cts <- cts[,-c(1,2)]

# Make the counts into a sparse matrix to save memory. 
cts <- as.matrix(cts) |> t()
cts <- as(cts, "CsparseMatrix")

# Note: you can still use base R data interaction on the sparse matrix. See:
cts[1:6, 1:6]
rownames(cts)

# Read in the metadata. 
meta <- read.csv("/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Grainger-Lab/Bioinformatics/Learning-Resources/Example-Data/Lung12_metadata_file.csv")

# Get the metadata into the correct format. 
rownames(meta) <- paste(meta$fov, meta$cell_ID, sep = "_")

# Make sure that the rownames of the metadata and the colnames of the counts match exactly. 
idx <- rownames(meta)
cts <- cts[,idx]

# Create the Seurat object. 
obj <- CreateSeuratObject(counts = cts, meta.data = meta, project = "training_obj")

# Skipping QC for now. 

# Normalize the data and find the most variable features. Apparently 2000-3000 is sufficient. 
obj <- NormalizeData(obj, normalization.method = "LogNormalize")
obj <- FindVariableFeatures(obj, method = "vst", nfeatures = 2000)

# Look for metadata variables that can predict the number of counts in a cell. 
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# If you see a variable that is positively correlated with nCount_RNA, then it 
# makes sense to regress it out. nFeature_RNA should be positively correlated 
# with nCount_RNA though (this makes intuitive sense), and should not be 
# regressed out, typically. We will regress it out here just for example's sake
# (so we can get used to working with the functions).

# Do some behind-the-scenes stuff to speed up the following steps. 
options(future.globals.maxSize = 2e9)
plan(strategy = "multisession", workers=4)

# Scale the variable features. Regress out the variable(s) you found above. 
obj <- ScaleData(obj, vars.to.regress = "nFeature_RNA")

# Run Principle Component Analysis (uses the scaled data) to do dimensional reduction.
obj <- RunPCA(obj, npcs = 50)

# Choose appropriate number of PCs. You want to choose a number past where the elbow is on this plot. 
ElbowPlot(obj, ndims = 50)

# Run the UMAP with the chosen PCs. 
obj <- RunUMAP(obj, dims = 1:20)

# Reverse those behind-the-scenes things. 
plan(strategy = "sequential")

# You now have the ability to view the UMAP. 
DimPlot(obj)

