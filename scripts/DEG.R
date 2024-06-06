#!/usr/bin/env Rscript

#.libPaths('/home/ash-maas/R/x86_64-pc-linux-gnu-library/4.4')

##import packages
require(limma)
require(RColorBrewer)
library(EnhancedVolcano)
library(biomaRt)


args <- commandArgs(trailingOnly=TRUE)

#setting input and output directories
#qc_results_dir<-'../output/pre-processing-output'
#inputfilepath<-file.path(qc_results_dir, "qc_data/")
#deg.out.dir<-"../output/deg-output/"
#outputplots.deg<-file.path(deg.out.dir, "plots")
#outputdata.deg<- file.path(deg.out.dir, "data")

#Import normalised expression and clinical data file
normalised_data<- read.delim(args[1],header=TRUE)

clinical.data<-read.delim(args[2],header=TRUE)


#1.5 Statistical analysis (DEG analysis)
#-------------------------------------------------------------------------------------------------------------------------------
#Create a design matrix.
design<-model.matrix(~0+factor(c(rep(c("HNL"), each=c(11)),rep(c("PSVD"), each = c(16))))+clinical.data$Sex)


#assign the column names.Here the cont_matrix represents in a matrix for the particular comparisons being made which indirectly
#represent the fold change for the particular comparison.
colnames(design)<-c("HNL","PSVD","gender")
cont_matrix<-makeContrasts(PSVDvsHNL=PSVD-HNL,
                           levels=design)

#Fit the expression matrix to a linear model.lmFit function of the "limma" package helps linearise the data with respect 
#to each gene and its corresponding fold change in the conditions i.e. row wise linearisation. This then allows for 
#selecting the genes with maximum fold change across conditions.
fit<-lmFit(normalised_data,design)

#Compute contrast
fit_contrast<-contrasts.fit(fit,cont_matrix)

#Bayes statistics of differential expression. Here the function eBayes performs t-Test and F-test to determine the p-value or significance
#by comparing means of the two conditions(t-test) and comapring two population variance(F-test).Additionally it also corrects for multiple
#testing by using q value(false detection rate) of 0.05 generally.
fit_contrast<-eBayes(fit_contrast)

#Generate a volcano plot to visualize the foldchanges across comparisons of conditions and identify the pattern of up- or down- regulation
#of the genes. Here we observe most of the genes getting down-regulated.
#volcanoplot(fit_contrast,xlab="log2Foldchange",ylab="-log10(pvalue)" )

#dataframe for foldchanges for psbd vs HNL comparison
psvdvshnl<-topTable(fit_contrast,coef = 'PSVDvsHNL',sort.by = "B",number=nrow(normalised_data),adjust.method ="BH")


#Summary of results(number of differentially expresses genes)
result<-data.frame(decideTests(fit_contrast))

#summary provides with number of genes which are either up-, down- or not significantly regulated across the comparisons.
summary(result)

##for converting entrez id to gene symbols 
psvdvshnl$entrez_id<-rownames(psvdvshnl)
mart <- useMart('ensembl',dataset="hsapiens_gene_ensembl")
filters <- listFilters(mart)
attributes <- listAttributes(mart)

gene.attributes <- getBM(mart=mart, values= psvdvshnl$entrez_id,
                         filters= "entrezgene_id",
                         attributes= c("entrezgene_id", "hgnc_symbol",
                                       "chromosome_name", "start_position",
                                       "end_position", "strand"))


deg.gene.data<-merge(gene.attributes, psvdvshnl, by.x ="entrezgene_id", by.y="entrez_id")

##keep unique gene symbol
deg.gene.data <- deg.gene.data[!duplicated(deg.gene.data$hgnc_symbol), ]

deg.gene.data.subset<-subset(deg.gene.data, abs(logFC) > 5.8 )


##volcano plot width = 20, height = 20
png('volcanoplot.png', width = 1000, height = 1000)
EnhancedVolcano(deg.gene.data, 
                lab = deg.gene.data$hgnc_symbol,
                x = 'logFC', y = 'adj.P.Val', 
                selectLab = deg.gene.data.subset$hgnc_symbol,
                pCutoff = 0.05, 
                FCcutoff = 1,
                ylab = bquote(~-log[10] ~ italic(adjusted.P)),
                xlab = bquote(log[2]~ foldchange),
                legendLabels = c("NS", expression(Log[2] ~ FC), "adjusted p-value", expression(adjusted ~ p - value ~ and
                                                                                               ~ log[2] ~ FC)),
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                ylim=c(0, 16),
                labSize = 7,
                axisLabSize = 25,
                captionLabSize = 25,
                legendLabSize = 25)
dev.off()


#export deg list with log fold change, p value, p adjusted value, B,..
write.table(psvdvshnl, "deg.data.txt", na ="", row.names=TRUE,  sep='\t', quote=FALSE)

