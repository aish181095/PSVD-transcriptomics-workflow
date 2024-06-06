#!/usr/bin/env Rscript

#.libPaths('/home/ash-maas/R/x86_64-pc-linux-gnu-library/4.4')

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
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)

##import module eigengene and module trait correlation data file from directory
MEs<-read.delim(args[1])
moduleTraitPvalue<-read.delim(args[2])

#subset module significantly correlating to the diagnosis of the paients 
moduleTraitpvalue_diagnosis<- subset(moduleTraitPvalue, Diagnosis <=0.05 )
MEs_diagnosis<-MEs[,rownames(moduleTraitpvalue_diagnosis)]


##calculatin correlation of module eigengenes to each other 
corr.mat<-cor(MEs_diagnosis, method = "pearson")

sampleTree = stats::hclust(dist(t(MEs_diagnosis)), method = "ward.D2");
plot(sampleTree, main = "Sample clustering to detect outliers", hang=-1,sub="", xlab="", cex.lab = 2, 
     cex.axis = 2, cex.main = 2)

#make node and edge table
node_corr_table<-data.frame(id = rownames(corr.mat))

edge_corr_table<-melt(corr.mat)
edge_corr_table <- edge_corr_table[-c(which(edge_corr_table$Var1 == edge_corr_table$Var2)),]
edge_corr_table<-subset(edge_corr_table, abs(value) > 0.5 ) ##cut -off for edge

##column in dataframe if edge represent positive or negative correlation.
edge_corr_table$edge_logical[edge_corr_table$value > 0] <- "positive"
edge_corr_table$edge_logical[edge_corr_table$value < 0] <- "negative"

edge_corr_table$value<- abs(edge_corr_table$value)

##check cytoscape connection
#cytoscapePing()
#if("cytargetlinker" %in% commandsHelp("")) print("Success: the CyTargetLinker app is installed") else print("Warning: CyTargetLinker app is not installed. Please install the CyTargetLinker app before proceeding.")

##export network to cytoscape
pathway_gene_nodes <- data.frame(id=node_corr_table$id,
                                 stringsAsFactors = FALSE)#node table


pathway_gene_edges <- data.frame(source=edge_corr_table$Var1,
                                 target=edge_corr_table$Var2,
                                 interaction=replicate(nrow(edge_corr_table), "contains"),
                                 weight = edge_corr_table$value,
                                 correlation = edge_corr_table$edge_logical,
                                 stringsAsFactors = FALSE)#edge table
#add column names to edge table
colnames(pathway_gene_edges)<-c("source", "target", "interaction", "weight", "correlation")

##export files
write.table(pathway_gene_nodes, "corrnet-node-table.txt", na ="", row.names=TRUE,  sep='\t', quote=FALSE)

write.table(pathway_gene_edges, "corrnet-edges-table.txt", na ="", row.names=TRUE,  sep='\t', quote=FALSE)

##create network in Cytoscape
#createNetworkFromDataFrames(pathway_gene_nodes, pathway_gene_edges, title="module", collection="wgcna")



