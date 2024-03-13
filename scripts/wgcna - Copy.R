##clear environment
remove(list = ls())# only run at the start of the script

if(!"knitr" %in% installed.packages()){
  install.packages("knitr")
}
library(knitr)
knitr:::input_dir()

##required libraries
library(AnnotationDbi)
library(preprocessCore)
library(GO.db)
library(readxl)
library(WGCNA)
library(rstudioapi)
library(dplyr)
library(biomaRt)
library(clusterProfiler)
library(flashClust)
library(pheatmap)
library(ComplexHeatmap)

# general config
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##set data and output directory
inputfilepath<-"../output/pre-processing-output/qc_data/"
outputfileplots<- '../output/wgcna-network-output/plots'#directory for output plots


#load normalised gene expression data
data <- read.delim(paste0(inputfilepath, 'normalised-data-afteroutlier.txt'),header=TRUE)
#transpose data
data<-as.data.frame(t(data))

# load clinical data 
clinical.data <- read.delim(paste0(inputfilepath, 'clinical-data-afteroutlier.txt'), header = T)


##check for minimal samples to run wgcna
if (nrow(data)|ncol(data) < 15){
  paste0("A minimum of 15 samples are required for WGCNA.")
}

#########################################
#Network construction 
#########################################

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.

enableWGCNAThreads(nThreads = 8)# Load the data saved in the first part
#lnames = load(file = "WGCNA-input.RData");
#The variable lnames contains the names of loaded variables.
#lnames


# Choose a set of soft-thresholding powers

powers = c(c(1:10), seq(from = 12, to=30, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(data,
                        powerVector = powers, 
                        verbose = 5, 
                        networkType = "signed")

# Plot scale free plot
pdf(paste0(outputfileplots, '/scale-free-plot.pdf'), width = 15, height = 10)
par(mar=c(5,6,4,1)+.1)
par(mfrow = c(1,2));
cex1 = 1.5;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab=expression(bold("Soft Threshold (power)")),
     cex.lab=1.2, cex.main=2,
     ylab=expression(bold("Scale Free Topology Model Fit, signed R^2")),type="n",
     main = expression(bold("Scale independence")));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], cex.lab=1.2,cex.main=2,
     xlab= expression(bold("Soft Threshold (power)")),ylab= expression(bold("Mean Connectivity")), type="n",
     main = expression(bold("Mean connectivity")));
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#soft threshold power selected
paste("Soft threshold power =",14)

##check the scale free topology of the network for the power selected
# When you have a lot of genes use the following code
k=softConnectivity(datE=data,
                   power=14, 
                   type="signed")

# Plot a histogram of k and a scale free topology plot
pdf(paste0(outputfileplots, '/scale-free-topology.pdf'), width = 15, height = 10)
par(mfrow=c(1,2))
hist(k, xlab = expression(bold("k")), ylab = expression(bold("Frequency")))
scaleFreePlot(k, main="Check scale free topology\n")
dev.off()


##Network construction
##Adjacency matrix calculation

softPower = 14;
adjacency = adjacency(data, 
                      power = softPower, 
                      type = "signed")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, 
                    TOMType = "signed");
dissTOM = 1-TOM

##save the TOM matrix
save(dissTOM, file = "WGCNA-input.RData")

lnames = load(file = "WGCNA-input.RData")

# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), 
                      method = "average");

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
plotDendroAndColors(geneTree,
                    dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

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
sizeGrWindow(7, 6)

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

####plot the new merged modules
pdf(paste0(outputfileplots, '/tree-dendrogram.pdf'), width = 15, height = 10)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

pdf(paste0(outputfileplots, '/module-eigengene.pdf'), width = 15, height = 10)
# Plot the relationships among the eigengenes 
plotEigengeneNetworks(MEs, "", marDendro = c(0, 4, 1, 2), 
                      marHeatmap = c(3, 4, 1, 2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()


moduleTraitCor = cor(MEs, clinical.data, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
#order the rows
moduleTraitCor<-moduleTraitCor[order(moduleTraitCor[,1], decreasing = TRUE),]
moduleTraitPvalue<-moduleTraitPvalue[rownames(moduleTraitCor),]

# Will display correlations and their p-values
textMatrix <- paste(formatC(moduleTraitCor, digits=2, format = "f"), "\n(",
                    formatC(moduleTraitPvalue, digits=1, format = "f"), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


#Heatmap for module-trait relationships 
pdf(paste0(outputfileplots, '/module-trait_relationship-all.pdf'), width = 40, height = 40)
# Will display correlations and their p-values
# Display the correlation values within a heatmap plot
par(mar=c(25,15,25,20)+.1)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinical.data),
               yLabels = rownames(moduleTraitPvalue),
               ySymbols = rownames(moduleTraitPvalue),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.5,
               cex.lab.x = 2,
               cex.lab.y = 2,
               cex.main=5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

moduleTraitPvalue<-as.data.frame(moduleTraitPvalue)
##plot for module correlating with diagnosis
subset_module_trait_relationship<-list()
for (names in colnames(moduleTraitPvalue)){
  subset_module_trait_relationship[[names]]<-subset(moduleTraitPvalue, moduleTraitPvalue[names] <=0.05)
}


##select desired clinical trait 
moduleTraitPvalue_diagnosis<-subset_module_trait_relationship[["Diagnosis"]]

moduleTraitCor_diagnosis<-moduleTraitCor[rownames(moduleTraitPvalue_diagnosis),]
##order the rows in descending order
moduleTraitCor_diagnosis<-moduleTraitCor_diagnosis[order(moduleTraitCor_diagnosis[,1], decreasing = TRUE),]

moduleTraitPvalue_diagnosis<-as.matrix(moduleTraitPvalue_diagnosis[rownames(moduleTraitCor_diagnosis),])

##Heatmap for modules significantly correlating to modules
diagnosis_matrix<-as.matrix(moduleTraitCor_diagnosis[,1])
colnames(diagnosis_matrix)<-"Diagnosis"
color.matrix<-as.matrix(1:15)
rownames(color.matrix)<-rownames(moduleTraitCor_diagnosis)

p1<- pheatmap(as.matrix(diagnosis_matrix[,1]),
         color = colorRampPalette(c("blue", "white", "red"))(50),
         labels_col = "Diagnosis",
         border_color = NA,
         fontsize_col = 9,
         angle_col = c("45"),
         legend = FALSE,
         number_color = "green",
         fontsize_number = 30,
         fontface = "bold",
         display_numbers = matrix(ifelse(moduleTraitPvalue_diagnosis[,1] < 0.05, "*", ""), nrow(moduleTraitPvalue_diagnosis)))
p2 <-  pheatmap(color.matrix,
           legend = FALSE,
           show_rownames = TRUE,
           border_color = NA,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           fontsize = 15,
           angle_col = "45",
           color= gsub('^..', '', rownames(moduleTraitCor_diagnosis)))
p3<- pheatmap(moduleTraitCor_diagnosis[,2:35],
         color = colorRampPalette(c("blue", "white", "red"))(50),
         name="correlation",
         border_color = NA,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         fontsize_col = 9,
         fontsize_row = 12,
         show_rownames = TRUE,
         angle_col = c("45"),
         fontsize_number = 30,
         fontface = "bold",
         number_color = "green",
         display_numbers = matrix(ifelse(moduleTraitPvalue_diagnosis[,2:35] < 0.05, "*", ""), nrow(moduleTraitPvalue_diagnosis)))

##save in the output directory 
pdf(paste0(outputfileplots, '/module-trait_relationship-subset.pdf'), width = 15, height = 8)
draw(p2+p1+p3)
dev.off()


# Define variable containing the diagnosis of the patients
diagnosis = as.data.frame(clinical.data$Diagnosis);
names(diagnosis) = "diagnosis"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(data, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(data, diagnosis, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(diagnosis), sep="");
names(GSPvalue) = paste("p.GS.", names(diagnosis), sep="");

# Create the starting data frame containing gene id , module they are present in, significance and correlation for each gene
geneInfo0 = data.frame(Gene.Symbol = colnames(data),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

##export module trait significance as a txt file
write.table(geneInfo0, '../output/wgcna-network-output/data/module-trait.txt', na ="", row.names=TRUE,  sep='\t', quote=FALSE)

##export module eigenegene as txt file 
write.table(MEs, "../output/wgcna-network-output/data/Module eigengene.txt", na ="", row.names=TRUE,  sep='\t', quote=FALSE)

##export module trait correlation significance matrix 
write.table(moduleTraitPvalue, "../output/wgcna-network-output/data/moduleTraitPvalue.txt", na ="", row.names=TRUE,  sep='\t', quote=FALSE)





