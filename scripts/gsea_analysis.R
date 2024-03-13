#!/usr/bin/env Rscript
 source("DEG.R")






##For the visualisation of gsea output open Cytoscape to view the network. 

##clear environment
remove(list = ls())# only run at the start of the script

if(!"knitr" %in% installed.packages()){
  install.packages("knitr")
}
library(knitr)
knitr:::input_dir()

##load required packages
library(clusterProfiler)
library(dplyr)
library(Hmisc)
library(msigdbr)
library(RCy3)
library(enrichplot)

# general config
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##set data and output directory
inputfilepath<-"../output/deg-output/data/"
outputfileplots<- '../../output/gsea-viz_output/plots'#directory for output plots

#import data: containing logFC, p-value, adjusted p-value and gene ids. 
psvdvshnl<- read.delim(paste0(inputfilepath, 'deg.data.txt'),header=TRUE, row.names = 1)


##create a genelist for GSEA
#ranking gene on p-value
psvdvshnl$ranking<-psvdvshnl$logFC*-log10(psvdvshnl$P.Value)
#ranking on product of logFC*-log10(adjusted p-value)
psvdvshnl$ranking_alt<-psvdvshnl$logFC*-log10(psvdvshnl$adj.P.Val)

#select the ranking: logFC*-log10(adjusted p-value)
geneList= psvdvshnl[,9]
#add entrex gene id as gene names
names(geneList) = as.character(rownames(psvdvshnl))
#sort on decreasing order
geneList = sort(geneList,decreasing = TRUE)


#import msigdb gene sets

kegg_gene_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
reactome_gene_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME")
wikipathways_gene_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:WIKIPATHWAYS")
bp_gene_sets<-msigdbr(species = "human", category = "C5", subcategory = "GO:BP")

#interchange the column name for bp_gene_sets "gs_description' and gs_name'
gs_name<- bp_gene_sets$gs_description
bp_gene_sets$gs_description<-NULL
names(bp_gene_sets)[names(bp_gene_sets)=="gs_name"]<-"gs_description"
bp_gene_sets$gs_name<-gs_name

##combine the WikiPathway, Reactome and KEGG gene sets
pathway_gene_sets<-rbind(bp_gene_sets, wikipathways_gene_sets, kegg_gene_sets, reactome_gene_sets)


##process the gene set as the input for the GSEA
pathway_gene_sets<- pathway_gene_sets %>%
  dplyr::select(gs_exact_source, entrez_gene, gs_description, gs_subcat) %>%
  dplyr::rename(ont = gs_exact_source, gene = entrez_gene, description = gs_description, source = gs_subcat )


##GSEA 
pathway_gsea<-GSEA(geneList,
                   TERM2GENE = pathway_gene_sets,
                   maxGSSize = 1000,
                   minGSSize = 30,
                   seed = TRUE,
                   by = "fgsea",
                   pAdjustMethod ="fdr",
                   pvalueCutoff = 0.05
                   
)
head(pathway_gsea,30)

#import the clusterprofiler enrichment output
pathway_names<-data.frame(pathway_gsea@result[["ID"]])


##create a list including the pathways enriched and their genesets
complete=list()
for (pathways in pathway_names$pathway_gsea.result...ID...){
  complete[[pathways]]<-data.frame(genes = pathway_gsea@geneSets[[pathways]],
                                   pathway_id = rep(c(pathways), each = nrow(complete[[pathways]])))
  
}

##Create a dataframe with pathways and genesets

pathway_data<-bind_rows(complete)

##calculate no of genes up and down regulated significantly: up : 1 and down:2 nd unchanged:0
psvdvshnl$names<-rownames(psvdvshnl)
psvdvshnl$condition[psvdvshnl$logFC >= 1 & psvdvshnl$adj.P.Val < 0.05] <- 1
psvdvshnl$condition[psvdvshnl$logFC <= -1 & psvdvshnl$adj.P.Val < 0.05] <- 2
psvdvshnl$condition[is.na(psvdvshnl$condition)]<- 0


##merging the two dataframes based on gene ids
merged_data<-merge(pathway_data, psvdvshnl, by.x = "genes", by.y = "names", all.x = TRUE)

#remove rows of NAs from the data
merged_data<-merged_data[complete.cases(merged_data),]

##create edge table
combined_pathways<-as.data.frame(merged_data[,1:2])

#counting frequencies
Counting<- data.frame(with(merged_data, unclass(table(pathway_id, condition))))

##renaming the columns
colnames(Counting)<-c("not-changed", "up-sig", "down-sig")
##calculating the total genes in a geneset
Counting$total<-Counting$`not-changed`+Counting$`up-sig`+Counting$`down-sig`
Counting$id<-rownames(Counting)

##Adding description for the pathway id
pathway_ont<-merge(Counting, pathway_gene_sets[,c(1,3,4)], by.x = "id", by.y = "ont", all.x = TRUE)
pathway_ont<-pathway_ont%>%
  distinct()

##Adding additonal columns
pathway_ont$logFC<-c(rep(c(NA), each = nrow(pathway_ont)))
pathway_ont$AveExpr<-c(rep(c(NA), each = nrow(pathway_ont)))
pathway_ont$t<-c(rep(c(NA), each = nrow(pathway_ont)))
pathway_ont$P.Value<-c(rep(c(NA), each = nrow(pathway_ont)))
pathway_ont$adj.P.Val<-c(rep(c(NA), each = nrow(pathway_ont)))
pathway_ont$B<-c(rep(c(NA), each = nrow(pathway_ont)))
pathway_ont$condition<-c(rep(c(NA), each = nrow(pathway_ont)))
pathway_ont$names<-c(rep(c(NA), each = nrow(pathway_ont)))
pathway_ont$entrez_id<-c(rep(c(NA), each = nrow(pathway_ont)))
rownames(pathway_ont)<-pathway_ont$id

##add NES score for pathway
NES_dataframe<-data.frame(id = pathway_gsea@result[["ID"]],
                          NES_score= pathway_gsea@result[["NES"]])
pathway_ont<-merge(pathway_ont, NES_dataframe, by.x="id", by.y="id", all.x=TRUE)

setsize_dataframe<-data.frame(id = pathway_gsea@result[["ID"]],
                              setSize= pathway_gsea@result[["setSize"]])
pathway_ont<-merge(pathway_ont, setsize_dataframe, by.x="id", by.y="id", all.x=TRUE)

##Adding additonal columns
psvdvshnl$id <- c(rep(c(NA), each = c(nrow(psvdvshnl))))
psvdvshnl$`up-sig` = c(rep(c(NA), each = c(nrow(psvdvshnl))))
psvdvshnl$`down-sig`<-c(rep(c(NA), each = c(nrow(psvdvshnl))))
psvdvshnl$`not-changed`<-c(rep(c(NA), each = c(nrow(psvdvshnl))))
psvdvshnl$total<-c(rep(c(NA), each = c(nrow(psvdvshnl))))
psvdvshnl$description<-c(rep(c(NA), each = c(nrow(psvdvshnl))))
psvdvshnl$source<-c(rep(c(NA), each = c(nrow(psvdvshnl))))
psvdvshnl$NES_score<-c(rep(c(NA), each = c(nrow(psvdvshnl))))
psvdvshnl$setSize<-c(rep(c(NA), each = c(nrow(psvdvshnl))))
psvdvshnl$ranking<-NULL
psvdvshnl$ranking_alt<-NULL
psvdvshnl$NES_score_condition<-c(rep(c(NA), each = c(nrow(psvdvshnl))))

#remove rows of NAs from the data
merged_data<-merged_data[complete.cases(merged_data),]

##create edge table
combined_pathways<-as.data.frame(merged_data[,1:2])

##node table
gene_id<-data.frame(id=unique(merged_data$genes))
gene_id$group<-c(replicate(nrow(gene_id), "Gene"))
pathway_id<-data.frame(id=pathway_ont$id,
                       group=pathway_ont$source)
node_table<-data.frame(id = rbind(gene_id, pathway_id))
colnames(node_table)<-c("id", "group")

#combine the pathway and gene table
cytoscape_table<-rbind(psvdvshnl, pathway_ont)

##calculating mutual nodes 
list_df<-split(combined_pathways, combined_pathways$pathway_id)
##split into multiple dataframes
list2env(list_df, envir= .GlobalEnv)

##mutual nodes


#pathway_ids<-pathway_names$pathway_gsea.result...ID...

pathway_enrich<-as_tibble(pathway_id$id)


rows= c(1:nrow(pathway_id))
times = nrow(pathway_id)
pathway_ids<-as.data.frame(pathway_id[rep(rows, times),])
pathway_ids$pathway_id_second<-pathway_enrich[rep(seq_len(nrow(pathway_enrich)), each=nrow(pathway_id)),]
pathway_ids<-as.data.frame(pathway_ids)

for (i in 1:nrow(pathway_ids)){
  
  
  pathway_ids$mutual_nodes[i]<-length(intersect(list_df[[pathway_ids[i,1]]][["genes"]], list_df[[unlist(pathway_ids[i,3])]][["genes"]]));
  
  pathway_ids$union_node[i]<-length(union(list_df[[pathway_ids[i,1]]][["genes"]], list_df[[unlist(pathway_ids[i,3])]][["genes"]]));
  
}

colnames(pathway_ids)<-c("source", "database", "target", "mutual_nodes", "union_nodes" )

setsize_table<-subset(pathway_ont, select = c("id", "setSize"))

setsize_match<-match(pathway_ids$source, setsize_table$id)
pathway_ids$source_setsize<-setsize_table[setsize_match,]

setsize_match_target<-match(pathway_ids$target$value, setsize_table$id)
pathway_ids$target_setsize<-setsize_table[setsize_match_target,]

pathway_ids$source_setsize$id<-NULL
pathway_ids$target_setsize$id<-NULL
pathway_ids$min_setsize<-pmin(pathway_ids$source_setsize, pathway_ids$target_setsize)

##overlap coefficient calculation 

pathway_ids$overlap_coef<- pathway_ids$mutual_nodes/pathway_ids$min_setsize$setSize


pathway_ids$jaccard_coef<- pathway_ids$mutual_nodes/pathway_ids$union_nodes
combined_constant = 0.5;

pathway_ids$combined_coef<- (combined_constant * pathway_ids$overlap_coef) + ((1-combined_constant) * pathway_ids$jaccard_coef)

#pathway_ids_subset<-subset(pathway_ids, overlap_coef > 0.5)

##make the edge table
pathway_ids_subset<-pathway_ids
pathway_ids_subset<-pathway_ids_subset[-c(which(pathway_ids_subset$source==pathway_ids_subset$target$value)),]
#pathway_ids_subset<-subset(pathway_ids_subset, overlap_coef$setSize == 1)

##adding an edge if the mutual nodes is above 50
pathway_ids$condition[pathway_ids$mutual_nodes> 50] <- 1
pathway_ids$condition[pathway_ids$mutual_nodes< 50] <- 0

pathway_ids<-subset(pathway_ids, condition ==1)

##nes score condition positive or negative
pathway_ont$NES_score_condition[pathway_ont$NES_score > 0] <- "positive"
pathway_ont$NES_score_condition[pathway_ont$NES_score < 0] <- "negative"

cytoscapePing()
if("cytargetlinker" %in% commandsHelp("")) print("Success: the CyTargetLinker app is installed") else print("Warning: CyTargetLinker app is not installed. Please install the CyTargetLinker app before proceeding.")

##Node table: Gene ID (Entrez) and Pathway Ids( KEGG, Reactome and WikiPathways)

pathway_gene_nodes <- data.frame(id=pathway_id$id,
                                 group= pathway_id$group,
                                 stringsAsFactors = FALSE)


# Edges : Gene-pathway associations
pathway_gene_edges <- data.frame(source=pathway_ids_subset$source,
                                 target=pathway_ids_subset$target$value,
                                 weight = pathway_ids_subset$overlap_coef,
                                 interaction=replicate(nrow(pathway_ids_subset), "contains"),
                                 stringsAsFactors = FALSE)

colnames(pathway_gene_edges)<-c("source", "target", "weight", "interaction" )
##create network in Cytoscape
createNetworkFromDataFrames(pathway_gene_nodes, pathway_gene_edges, title="GSEA pathway nodes cluster-jaccard", collection="GSEA")


cytoscape_table$description<-tolower(cytoscape_table$description)
cytoscape_table$description<-capitalize(cytoscape_table$description)

pathway_ont$description_new<-gsub("^.*?_","_", pathway_ont$description)
pathway_ont$description_new<-gsub("_", " ", pathway_ont$description_new)
pathway_ont$description_new<-tolower(pathway_ont$description_new)
pathway_ont$description_new<-capitalize(pathway_ont$description_new)
#pathway_ont$NES_score<-abs(pathway_ont$NES_score)


##load node table data 
#loadTableData(cytoscape_table)
rownames(pathway_ont)<-pathway_ont$id
loadTableData(pathway_ont)
#write.table(NES_dataframe, "nes_data.txt", na ="", row.names=FALSE,  sep='\t', quote=FALSE)

rownames(pathway_id)<-pathway_id$id
loadTableData(pathway_id)

##set the node shape 
getNodeShapes()   # diamond, ellipse, trapezoid, triangle, etc.
column <- 'group'
values <- c('CP:KEGG', 'CP:REACTOME', 'CP:WIKIPATHWAYS', 'GO:BP')
shapes <- c( 'ELLIPSE' , 'ELLIPSE', 'ELLIPSE', 'ELLIPSE')
setNodeShapeMapping(column, values, shapes)

##Adding the pie chart to pathway nodes to visualise the down-regulated, not -changed and the up-regulated genes.
setNodeCustomPieChart(c("down.sig","not.changed","up.sig"), colors = c('#334CFF', '#FFFFFF', '#FD5903') )

lockNodeDimensions(TRUE)

##set the node size
size<-c(100,100,100,100)
setNodeSizeMapping('group', mapping.type = "d", values, size)

#set the node label: only set for the pathway nodes
setNodeLabelMapping("description_new")

#set the node position
setNodeCustomPosition(nodeAnchor = "C", graphicAnchor = "C", justification = "c")

#set the node label font size
setNodeFontSizeDefault(20)

#set the font style of the node label
setNodeFontFaceDefault("Arial,Rounded,Bold,45")

#set the node border width
border_width=c(20,20,20,20, 20)
setNodeBorderWidthMapping('group', mapping.type = "d", values, border_width)

#set the node border color: for visualising the pathway databases
border_color=c("#FD5903", "#0026D3" )
setNodeBorderColorMapping('NES_score_condition', mapping.type = "d", values, border_color)

#set the edge color
setEdgeColorDefault("#D7DBDD")

deleteDuplicateEdges()
