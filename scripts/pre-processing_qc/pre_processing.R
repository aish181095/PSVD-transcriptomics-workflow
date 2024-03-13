#!/usr/bin/env Rscript

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

# general config

pre_processing<-function(raw_data, annot){

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))##setting working directory path

#bgxfile <- 'Annot/HumanHT-12_V4_0_R2_15002873_B.bgx'##Illumina beadchip 
targetsfile <- 'targets.txt'
#outputfileplots<- '../../output/pre-processing-output/plots'#directory for output plots
#outputfiledata<-'../../output/pre-processing-output/qc-data'#directory for qc data 
options(scipen = 99)

# read in the data and convert the data to an EListRaw object, which is a data object for single channel data
x <- read.table(raw_data,
                header = TRUE, sep = '\t', stringsAsFactors = FALSE)
newnames <- colnames(read.table(raw_data,
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE))
detectionpvalues <- x[,grep('Detection.Pval', colnames(x))]
x <- x[,-grep('Detection.Pval', colnames(x))]
probes <- x$ID_REF
x <- data.matrix(x[,2:ncol(x)])
rownames(x) <- probes

#subset only HNL and PSVD
x<-x[,1:32]
detectionpvalues<-detectionpvalues[,1:32]

##add rownames
HNL_names<-c(sprintf("HNL%02d", seq(1,14), "label"))
INCPH_names<-c(sprintf("PSVD%02d", seq(1,18), "label"))

colnames(x)<-c(HNL_names,INCPH_names)


# read in annotation
annot <- illuminaio::readBGX(annot)$probes
annot <- annot[,which(colnames(annot) %in% c('Source','Symbol','Transcript','ILMN_Gene','RefSeq_ID',
                                             'Entrez_Gene_ID','Symbol','Protein_Product','Probe_Id','Probe_Type',
                                             'Probe_Start','Chromosome','Probe_Chr_Orientation','Probe_Coordinates',
                                             'Cytoband', 'Definition', 'Ontology_Component', 'Ontology_Process',
                                             'Ontology_Function', 'Synonyms'))]
annot <- annot[which(annot$Probe_Id %in% rownames(x)),]
annot <- annot[match(rownames(x), annot$Probe_Id),]

##make target file
targetinfo<-data.frame(IDATfile= colnames(x),
                       Group = c(rep(c("Healthy"), each = c(14)),rep(c("PSVD"), each = c(18))))

# update the target info
rownames(targetinfo) <- targetinfo$IDATfile
x <- x[,match(rownames(targetinfo), colnames(x))]
if (!all(colnames(x) == rownames(targetinfo)))
  stop('Target info is not aligned to expression data - they must be in the same order')

# create a custom EListRaw object
project <- new('EListRaw')
project@.Data[[1]] <- 'illumina'
project@.Data[[2]] <- targetinfo
project@.Data[[3]] <- annot
project@.Data[[3]] <- NULL
project@.Data[[4]] <- x
project@.Data[[5]] <- NULL
project$E <- x
project$targets <- targetinfo
project$genes <- annot
project$genes <- NULL
project$other$Detection <- detectionpvalues


# for BeadArrays, background correction and normalisation are handled by a single function: neqc()
# this is the same as per Agilent single colour arrays
#
# perform background correction on the fluorescent intensities
#   'normexp' is beneficial in that it doesn't result in negative values, meaning no data is lost
#   for ideal offset, see Gordon Smyth's answer, here: https://stat.ethz.ch/pipermail/bioconductor/2006-April/012554.html
# normalize the data with the 'quantile' method, to be consistent with RMA for Affymetrix arrays

project.bgcorrect.norm <- neqc(project, offset = 16)

# filter out control probes, those with no symbol, and those that failed
annot <- annot[which(annot$Probe_Id %in% rownames(project.bgcorrect.norm)),]
project.bgcorrect.norm <- project.bgcorrect.norm[which(rownames(project.bgcorrect.norm) %in% annot$Probe_Id),]
annot <- annot[match(rownames(project.bgcorrect.norm), annot$Probe_Id),]
project.bgcorrect.norm@.Data[[3]] <- annot
project.bgcorrect.norm$genes <- annot
Control <- project.bgcorrect.norm$genes$Source=="ILMN_Controls"
NoSymbol <- project.bgcorrect.norm$genes$Symbol == ""
isexpr <- rowSums(project.bgcorrect.norm$other$Detection <= 0.01) >= 5
project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]
dim(project.bgcorrect.norm)
dim(project.bgcorrect.norm.filt)


# remove annotation columns we no longer need
project.bgcorrect.norm.filt$genes <- project.bgcorrect.norm.filt$genes[,c(
  'Probe_Id', 'Entrez_Gene_ID',
  'Definition','Ontology_Component','Ontology_Process','Ontology_Function',
  'Chromosome','Probe_Coordinates','Cytoband','Probe_Chr_Orientation',
  'RefSeq_ID','Entrez_Gene_ID','Symbol')]
head(project.bgcorrect.norm.filt$genes)


# remove annotation columns we no longer need
project.bgcorrect.norm.filt$genes <- project.bgcorrect.norm.filt$genes[,c(
  'Definition','Ontology_Component','Ontology_Process','Ontology_Function',
  'Cytoband','Probe_Chr_Orientation',
  'RefSeq_ID','Symbol', 'Entrez_Gene_ID')]
head(project.bgcorrect.norm.filt$genes)

# summarise across genes by mean
# ID is used to identify the replicates
project.bgcorrect.norm.filt.mean <- avereps(project.bgcorrect.norm.filt,
                                            ID = project.bgcorrect.norm.filt$genes$Entrez_Gene_ID)
dim(project.bgcorrect.norm.filt.mean)

##expression data after normalisation
expression_data<-as.data.frame(project.bgcorrect.norm.filt.mean[["E"]])

##Detection probe
detection_probe<-as.data.frame(project.bgcorrect.norm.filt.mean[["other"]][["Detection"]])


#############
#Plots
#############

#plot the expression data (before normalisation)
width= 20;
height=15;
pdf(paste0(outputfileplots, '/before-normalisation.pdf'), width = width, height = height)
box_plot_data<-as.data.frame(log2(x))
par(mar=c(5,6,4,1)+.1)
colors=c(rep(2,14),rep(4,18))
boxplot(box_plot_data, col= colors,main="Before background correction and normalisation",
        xlab="Samples",ylab=expression('log'[2]*'(expression)'),cex=0.5, par(cex.lab=1.5),par(cex.axis=1.5))
colors_names<-c(2,4)
legend("topright",c("Healthy","PSVD"),fill=colors_names, cex=1.5)
dev.off()

#plot the detection score (p values indicating how significant each read is)
pdf(paste0(outputfileplots, '/before-normalisation-probe.pdf'), width = width, height = height)
boxplot(detectionpvalues,main="Detection score (pvalue) before background correction ", col=colors,
        xlab="Samples",ylab="P value",cex=0.8)
dev.off()


#plot the expression data (before normalisation)
pdf(paste0(outputfileplots, '/after-normalisation.pdf'), width = width, height = height)
par(mar=c(5,6,4,1)+.1)
colors=c(rep(2,14),rep(4,18))
boxplot(expression_data, col= colors,main="After background correction and normalisation",
        xlab="Samples",ylab=expression('log'[2]*'(expression)'),cex=0.5, par(cex.lab=1.5),par(cex.axis=1.5))
colors_names<-c(2,4)
legend("topright",c("Healthy","PSVD"),fill=colors_names, cex=1.5)
dev.off()

#plot the detection score (p values indicating how significant each read is)
pdf(paste0(outputfileplots, '/after-normalisation-probe.pdf'), width = width, height = height)
boxplot(detection_probe,main="Detection score (pvalue) after background correction ", col=colors,
        xlab="Samples",ylab="P value",cex=0.8)
dev.off()


##export data to output directory
write.table(expression_data, "normalised.data.txt", na ="", row.names=TRUE,  sep='\t', quote=FALSE)
}

