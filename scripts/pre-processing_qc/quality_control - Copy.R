if(!"knitr" %in% installed.packages()){
  install.packages("knitr")
}
library(knitr)
knitr:::input_dir()

##import packages
require(limma)
require(RColorBrewer)
require(PCAtools)
library(illuminaio)
library(WGCNA)
library(readxl)
library(dendextend)

# general config
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##set data and outpuit directory
inputfilepath<-"../../output/pre-processing-output/qc_data/"
outputfileplots<- '../../output/pre-processing-output/plots'#directory for output plots
outputfiledata<-'../../output/pre-processing-output/qc-data/'#directory for qc data 

# Load the data
data <- read.delim(paste0(inputfilepath, 'normalised.data.txt'),header=TRUE)

# remove misdiagnosed patient (PSVD17)
data <- subset(data, select = -(PSVD17))

# load clinical data and keep only relevant columns
clinical.data <- read_excel("../../Data/clinical-data/clinical_data_intial.xlsx")
clinical.data <- subset(clinical.data, select = -c(id, NHC, `alt id`, Ageatdx, `Ascites at diagnosis` , IPH, CH, HV, `Which associated disease`))


# remove CIRR samples
clinical.data.filt <- as.data.frame(clinical.data[clinical.data$Diagnosis != "CIRR",])
rownames(clinical.data.filt) <- colnames(data)

# remove data with missing sex information (HNL07, HNL13, HNL14)
remove <- c("HNL07","HNL13","HNL14")
clinical.data.filt <- clinical.data.filt[!(row.names(clinical.data.filt) %in% remove),] 
data.filt <- data[,!(colnames(data) %in% remove)]

#convert categorical variables to numeric levels 
#first convert the columns to factors then to numeric values
clinical.data.filt$Diagnosis <- as.numeric(unclass(as.factor(clinical.data.filt$Diagnosis)))
clinical.data.filt$Sex <- as.numeric(unclass(as.factor(clinical.data.filt$Sex)))
clinical.data.filt$`Associated disease` <- as.numeric(unclass(as.factor(clinical.data.filt$`Associated disease`)))
#clinical.data.filt$`Ascites at diagnosis` <- as.numeric(unclass(as.factor(clinical.data.filt$`Ascites at diagnosis`)))
clinical.data.filt$HVPG <- as.numeric(clinical.data.filt$HVPG)
clinical.data.filt$WHVP <- as.numeric(clinical.data.filt$WHVP)
#clinical.data.filt$`Platelet count` <- as.numeric(clinical.data.filt$`Platelet count`)
clinical.data.filt$`Total bilirubin` <- as.numeric(clinical.data.filt$`Total bilirubin`)
clinical.data.filt$`Direct bilirubin at diagnosis` <- as.numeric(clinical.data.filt$`Direct bilirubin at diagnosis`)
clinical.data.filt$`Spleen size` <- as.numeric(clinical.data.filt$`Spleen size`)
clinical.data.filt$`Liver elastography` <- as.numeric(clinical.data.filt$`Liver elastography`)
rownames(clinical.data.filt) <- colnames(data.filt)


# Data quality check 
gsg = goodSamplesGenes(data, verbose = 3);
gsg$allOK

# remove rows with all low expression
data.filt <- data.filt[rowSums(data.filt) > 0.5, , drop=TRUE]
data.filt <- as.data.frame(t(as.matrix(data.filt)))

##Sample clustering visualisation for outliers
sampleTree = stats::hclust(dist(data.filt), method = "ward.D2");
groupCodes <- c(rep("Healthy",11), rep("PSVD",4), rep("OUTLIER",1), rep("PSVD",12))
colorCodes <- c(Healthy="black", PSVD="blue", OUTLIER="red")
dend <- as.dendrogram(sampleTree, main = "Sample clustering to detect outliers", hang=-1,sub="", xlab="", cex.lab = 2, cex.axis = 2, cex.main = 2)

# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]

# Plotting sample clustering
pdf(paste0(outputfileplots, '/sample-clustering.pdf'), width = 20, height = 10)
par(mar=c(5,6,4,1)+.1)
plot(dend,main = "Sample clustering to detect outliers",sub="", ylab="Height", cex.lab = 2, cex.axis = 2, cex.main = 2)+
scale_y_discrete(expand = c(0, 1))
dev.off()

#remove outlier sample (PSVD05)
remove <- c("PSVD05")
clinical.data.filt <- clinical.data.filt[!(row.names(clinical.data.filt) %in% remove),] 
data.filt <- data.filt[!(rownames(data.filt) %in% remove),]


##Sample clustering visualisation for outliers
sampleTree = stats::hclust(dist(data.filt), method = "ward.D2");
groupCodes <- c(rep("Healthy",11), rep("PSVD",16) )
colorCodes <- c(Healthy="black", PSVD="blue")
dend <- as.dendrogram(sampleTree, main = "Sample clustering to detect outliers", hang=-1,sub="", xlab="", cex.lab = 2, cex.axis = 2, cex.main = 2)


# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]

# Plotting the new dendrogram
pdf(paste0(outputfileplots, '/sample-clustering-after outlier removed.pdf'), width = 20, height = 10)
par(mar=c(5,6,4,1)+.1)
plot(dend,main = "Sample clustering after removing outliers",sub="", ylab="Height", cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()


##Sample and clinical data clustering

# Re-cluster samples
sampleTree2 = stats::hclust(dist(data.filt), method = "ward.D2")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(clinical.data.filt, signed = FALSE);

#Sample clustering and clinical data visualisation

pdf(paste0(outputfileplots, '/sample dendrogram and train heatmap.pdf'), width = 15, height = 15)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(clinical.data.filt), 
                    main = "Sample dendrogram and trait heatmap",
                    marAll = c(1, 25, 3, 1),
                    cex.colorLabels=1.0,
                    cex.dendroLabels = 1.5,
                    addGuide = FALSE)
dev.off()

##export data without outliers
write.csv(as.data.frame(t(data.filt)),file.path(outputfiledata, "norm-data-final.csv"))

write.table(as.data.frame(t(data.filt)), "C:\\Users\\P70069388\\Documents\\PVH data\\Illumina-processing\\git-repositories\\PSVD-Workflow\\output\\pre-processing-output\\qc_data\\normalised-data-afteroutlier.txt", na ="", row.names=TRUE,  sep='\t', quote=FALSE)
write.table(clinical.data.filt, "C:\\Users\\P70069388\\Documents\\PVH data\\Illumina-processing\\git-repositories\\PSVD-Workflow\\output\\pre-processing-output\\qc_data\\clinical-data-afteroutlier.txt", na ="", row.names=TRUE,  sep='\t', quote=FALSE)
