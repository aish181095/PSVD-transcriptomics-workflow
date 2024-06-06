cytoscapePing()
if("cytargetlinker" %in% commandsHelp("")) print("Success: the CyTargetLinker app is installed") else print("Warning: CyTargetLinker app is not installed. Please install the CyTargetLinker app before proceeding.")


##import data
##Node table: Gene ID (Entrez) and Pathway Ids( KEGG, Reactome and WikiPathways)
pathway_gene_nodes <-  read.delim("cytoscape-node-table.txt", header=TRUE, row.names = 1)


# Edges : Gene-pathway associations
pathway_gene_edges <- read.delim("cytoscape-edges-table.txt", header=TRUE, row.names = 1)

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

##save cytoscape session
session.file <- "gsea_visualisation.cys"
saveSession(session.file)