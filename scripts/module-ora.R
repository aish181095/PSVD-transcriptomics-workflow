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
library(msigdbr)
library(enrichplot)
library(Hmisc)
library(ggplot2)

library(org.Hs.eg.db)

args <- commandArgs(trailingOnly=TRUE)

##import module trait significance data
geneInfo0<- read.delim(args[1],header=TRUE)

geneInfo0$Gene.Symbol<-as.character(geneInfo0$Gene.Symbol)

##creating all module color dataframes
list_df <- split(geneInfo0, geneInfo0$moduleColor) #split the dataset into a list of datasets based on the value of iris$Species
list2env(list_df, envir= .GlobalEnv) #split the list into separate datasets

#msigdb gene sets
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

##create input (term2name and term2gene dataframe for enricher function)
pathway2gene<-pathway_gene_sets[,c("ont", "gene")]
pathway2name<-pathway_gene_sets[,c("ont", "description")]

##data wrangling for pathway term
pathway2name$description<-gsub("^.*?_","_", pathway2name$description)#removes "GO_" before the term
pathway2name$description<-gsub("_", " ", pathway2name$description)#removes "_" between terms
pathway2name$description<-gsub("^\\s+|\\s+$","",pathway2name$description)#removes blank space before the term
pathway2name$description<-tolower(pathway2name$description)#converts term to small case
pathway2name$description<-capitalize(pathway2name$description)#converts first letter of term to capital

##ora for each modules 
for (samples in c("grey60", "red", "black","turquoise", "yellow", "darkturquoise","lightyellow", "blue", "brown")) {
  ora <- enricher(
    gene = get(samples)$Gene.Symbol, # A vector of your genes of interest
    pvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,# Can choose a FDR cutoff
    pAdjustMethod = "fdr", # Method to be used for multiple testing correction
    # A vector containing your background set genes
    # The pathway information should be a data frame with a term name or
    # identifier and the gene identifiers
    TERM2GENE = pathway2gene,
    TERM2NAME = pathway2name,
    qvalueCutoff = 0.5,
    universe = geneInfo0$Gene.Symbol
  )
  newname<- paste(samples,"ora", sep ="_")
  assign(newname, ora)
  
  
}

for (samples in c("green", "salmon", "cyan", "lightcyan", "tan")) {
  ora <- enrichGO(
    gene = get(samples)$Gene.Symbol, # A vector of your genes of interest
    pvalueCutoff = 0.01,
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENTREZID',
    minGSSize = 10,
    maxGSSize = 500,# Can choose a FDR cutoff
    pAdjustMethod = "none", # Method to be used for multiple testing correction
    qvalueCutoff = 0.5,
    universe = geneInfo0$Gene.Symbol
  )
  newname<- paste(samples,"ora", sep ="_")
  assign(newname, ora)
}


for (samples in c("grey60_ora", "red_ora", "black_ora", "darkturquoise_ora","blue_ora", "brown_ora", "green_ora", "lightcyan_ora", "salmon_ora", "cyan_ora", "yellow_ora", "turquoise_ora", "lightyellow_ora")){
  if (length(get(samples)$ID)== 0){
    paste0("No terms enriched.")
  }
  else if (length(get(samples)$ID)<=4) {
    #dot plot 
    pdf(file=paste(samples, "_dotplot.pdf", sep=""), paper="a4")
    print(dotplot(get(samples), showCategory=45))
    dev.off()
  }
  else if (length(get(samples)$ID)>4)  {
    ##pairwise term calculation
    pairwise.calc<-pairwise_termsim(get(samples))
    newname.1<- paste(samples,"pairwise", sep ="_")
    assign(newname.1, pairwise.calc)
    
    ##make and save tree plot
    pdf(file=paste(gsub("_[^_]+$", "",newname.1), "_treeplot.pdf", sep=""), width = 25, height = 8)
    print(treeplot(get(newname.1), showCategory=45, cluster.params = list(method = "ward.D2", label_words_n = 0), label_format = 5, fontsize =6)+ 
            theme(legend.position="left"))
    dev.off()
  }
}