##required libraries
library(AnnotationDbi)
library(preprocessCore)
library(GO.db)
library(readxl)
library(WGCNA)
library(dplyr)
library(biomaRt)
library(clusterProfiler)
library(flashClust)
library(pheatmap)
library(ComplexHeatmap)


lnames = load(file = "WGCNA-input.RData")


# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), 
                      method = "average");

#load normalised gene expression data
data <- read.delim('/home/ash-maas/Documents/PSVD-transcriptomics-workflow/output/pre-processing-output/qc_data/normalised-data-afteroutlier.txt' , header=TRUE)
#transpose data
data<-as.data.frame(t(data))

# load clinical data 
clinical.data <- read.delim('/home/ash-maas/Documents/PSVD-transcriptomics-workflow/output/pre-processing-output/qc_data/clinical-data-afteroutlier.txt', header = T)
# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.07)

#set minimum module size
minModuleSize = 150;


# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, 
                            distM=dissTOM, method="hybrid", 
                            minClusterSize = minModuleSize, 
                            deepSplit = 4);
table(dynamicMods)

##assigning colors to modules
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

##Visualising the dendogram of the network and the modules detected
#plotDendroAndColors(geneTree,
#                    dynamicColors, "Dynamic Tree Cut", 
#                    dendroLabels = FALSE, 
#                    hang = 0.03, 
#                    addGuide = TRUE, 
#                    guideHang = 0.05, 
#                    main = "Gene dendrogram and module colors")

##cluster modules together that their expression is similar
MEList = moduleEigengenes(data, colors = dynamicColors)
# Calculate dissimilarity of module eigengenes
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
MECorr = cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")


##plot clustering of modules
#sizeGrWindow(7, 6)

plot(METree, 
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")
abline(h=0, col = "red")



#merging close modules
merge <- WGCNA :: mergeCloseModules(data, dynamicColors, cutHeight = 0,
                                    verbose = 3) 

# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs

nGenes=ncol(data)
nSamples=nrow(data)
MEs0=moduleEigengenes(data,moduleColors)$eigengenes
MEs=orderMEs(MEs0)