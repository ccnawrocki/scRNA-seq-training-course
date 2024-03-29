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

# Work through these steps on your own. This will take you a while and you will 
# struggle. You do not have to finish it today. 

# Start with reading the data. HINT: use fread() function from the data.table
# package. 

# Make sure all of the cell IDs match in the metadata and in the counts table.

# Next, do QC for each sample individually. 

# Combine the samples. HINT: use the merge() function.

# Filter out genes that do not appear in at least 10 cells. HINT: use rowSums()
# function. 

# Seurat standard workflow: Normalizing and scaling, clustering, etc. 

# Save the data as .Rds and as .h5Seurat

# Differential expression analysis to find cluster markers. 

# Differential expression analysis with pseudobulking to compare clusters 10 and 13. 
