### Training Day 5 ###
## By Cole Nawrocki ##

## Topic: Differential expression analysis ## 

# What is it DE analysis? 
# We want to identify genes that are expressed at significantly different levels 
# across two states. Let's say that we are interested in pax6 and we have some 
# scRNA-seq data with cells from a control group and an experimental group. We 
# can take all the cells in a control group and call that a sample. We can take 
# all the cells in the experimental group and call that a sample. Then, we can 
# find the mean expression level of pax6 in each sample and run a t test to 
# determine whether pax6 is expressed differently across the samples. If pax6 is 
# significantly different across the samples, then this would provide evidence 
# that the experimental treatment had some effect on the cells' expression of 
# the gene. However, how do you choose genes to look into? With computing power, 
# do do not have to choose. Just do the test for every gene and find the ones
# that are significant. This process is called differential expression analysis.

# Some more info: 
# When you do DEA, you do not typically use a t test, since that test has 
# specific assumptions that must be met. Often you use the wilcoxon rank-sum
# test, since it is a non-parametric test and therefore makes no assumptions 
# about the distribution of the population of cells. Furthermore, while taking
# cells individually as observations works, it is often much better to group 
# cells from the same sample together and pool their counts. This practice is 
# called "pseudobulking." Before scRNA-seq, there was (and still is) bulk 
# RNA-seq. Here, we are turning the scRNA-seq into bulk RNA-seq (sort of), so we
# call its pseudo-bulk.  

# Let's work through an example. We will use the briggs stage 16 data. 

# Load packages
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(pheatmap)

# Read in data
st16 <- LoadSeuratRds("/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Grainger-Lab/Bioinformatics/Data/Briggs-Data/Rds-files/briggs_data_stage16.Rds")
DimPlot(st16, label = T) + NoLegend()

# What makes cluster 3 unique? We can use DEA to find out. All cells in cluster 
# 3 will be a group and all the other cells will be another group. 

# Single-cell DEA method:
# The easiest way to do this is with Seurat's FindMarker's function. 
cl3_markers <- FindMarkers(st16, ident.1 = 3, logfc.threshold = 1, only.pos = T, test.use = "wilcox")

# The function knows that I am talking about cluster 3, because I used the 
# SetIdents() function earlier. It is very important to either reset the default
# identifiers with the SetIdents() function to the correct meta data column OR
# explicitly include the group.by argument with the correct meta data column if 
# you are going to compare groups from a different meta data column. For 
# instance, if we were going to compare groups in the "RNA_snn_res.0.1" column, we would
# do one of the following:

# Comparing clusters 0 and 1 in this column. 
# SetIdents(st16, value = "RNA_snn_res.0.1")
# test <- FindMarkers(st16, ident.1 = 0, ident.2 = 1, test.use = "wilcox", logfc.threshold = 1)
# OR
# test <- FindMarkers(st16, group.by = "RNA_snn_res.0.1", ident.1 = 0, ident.2 = 1, test.use = "wilcox", logfc.threshold = 1)

# After looking at the table, it is clear that tfap2b is a strong marker for 
# cluster 3, compared to the others. 

# Visualize individual genes like this:
FeaturePlot(st16, "tfap2b", order = T)

# It is easy to do this for every cluster: 
all_markers <- FindAllMarkers(st16, logfc.threshold = 1, only.pos = T, test.use = "wilcox")

# We can summarize marker genes for each cluster like this: 
top_10 <- all_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC)
DoHeatmap(st16, features = top_10$gene)
# You will want to make the plot prettier and save the markers in a csv. 

# You may be wondering: what is LFC? It is log fold change. This is the average 
# of group 1 divided by the average of group 2, plugged into the log2 function. 
# This means that when log2FC is very positive, the average for group 1 is 
# higher. When the log2FC is very negative, the average for group 2 is higher. 
# Above, I only wanted genes that were strong, positive markers for group 1, 
# which was each cluster. If you were doing DEA across experimental groups, 
# you would want negative markers too. 

# You may also be wondering: what is p_val_adj? It is p-value adjusted. 
# Essentially, when you do this many tests, there is an inevitable probability 
# of type I error. We essentially add this probability to each p-value in order 
# to make sure that none of our significant genes are just showing up due to 
# that error. The sum of this probability and the p-value is the p-value 
# adjusted and is de facto. 

# Pseudobulk DEA method: DESeq2
# Let's focus in on cluster 3 again. 
agg_counts <- AggregateExpression(st16, group.by = c("louvain_20_clusters", "Library_name"))[["RNA"]]

# Now the data looks like this. 
agg_counts[1:4,1:4]

# Let's make metadata
agg_meta <- data.frame(row.names = colnames(agg_counts))
agg_meta$cluster <- str_split(rownames(agg_meta), pattern = "_", simplify = T)[,1]
agg_meta$cluster3 <- ifelse(agg_meta$cluster == "g3", yes = "cluster3", no = "other")

# Now, we can use the DESeq2 model. This is a fancy RNA sequencing ata model 
# that will work better than the wilcoxon rank sum test. 
dds <- DESeqDataSetFromMatrix(countData = agg_counts, 
                              colData = agg_meta, 
                              design = ~cluster3)

dds <- DESeq(dds)
res <- results(dds) |> as.data.frame() |> na.omit()

# Visualize with a volcano plot. It looks like a volcano erupting. The genes in
# the top corners are "interesting."
res$delabel<-'Not Enriched'
res$delabel[res$padj<0.05 & res$log2FoldChange > 1]<-"De-Enriched"
res$delabel[res$padj<0.05 & res$log2FoldChange < -1]<-"Enriched"
res <- res[order(res$delabel),]
res$lbl<-NA
res$lbl[res$delabel != 'Not Enriched'] <- rownames(res)[res$delabel != 'Not Enriched']
v.plot<-ggplot(data=res,aes(x=log2FoldChange, y=-log10(padj),col=delabel, label=lbl)) +
  geom_point() + theme_minimal() + geom_text_repel(max.overlaps = 50) +
  labs(title="Enriched Genes in Cluster 3") + 
  scale_color_manual(values=c("De-Enriched"="red", "Enriched"="blue","Not Enriched"="gray")) +
  geom_vline(xintercept=c(-1,1),col="black") + 
  geom_hline(yintercept=-log10(0.05),col="black")
print(v.plot)

# Some other ways to visualize:
# DESeq2 plotCounts function:
plotCounts(dds, "tfap2b", intgroup = "cluster3")
plotCounts(dds, "sox10", intgroup = "cluster3")
plotCounts(dds, "krt", intgroup = "cluster3")

# Heatmap:
genes_oi <- res[(res$padj < 0.05) & (abs(res$log2FoldChange) > 4) & order(res$log2FoldChange),] |> rownames()
norm <- counts(dds, normalized = TRUE)
sig_counts <- norm[genes_oi, ]

pheatmap(sig_counts, 
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         annotation = agg_meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 5, 
         height = 20,
         fontsize_col = 3)  

# Final notes: 
# There are many other popular models/methods. These include edgeR, MAST, 
# and glmGamPoi. 
# The example above will not bridge you through every problem that you encounter
# by any means. You will have to get creative and READ papers/vignettes about 
# how to do certain comparisons with certain methods. This example is to get 
# you started with the basics and to get you thinking about how this works.
# The pseudobulking tends to be more valuable for comparing two large groups as 
# opposed to identifying markers for clusters. For example, the pseudobulking 
# would likely be my choice of method for comparing stage 18 to stage 22 or 
# mutant to WT. 

