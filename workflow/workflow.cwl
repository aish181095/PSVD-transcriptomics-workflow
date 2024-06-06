cwlVersion: v1.0
class: Workflow

inputs:  
 raw_counts:
  type: File
  label: Microarray data
  doc: Illumina beadchip raw gene expression data -GSE77627
  format: http://edamontology.org/data_2603 # Expression data
   
 bgxfile:
  type: File
  label: Illumina beadchip array file
  doc: Illumina HumanHT-12 WG-DASL V4.0 R2 expression beadchip 
  format: http://edamontology.org/data_2526 #text data
    
 clinical_data:
  type: File
  label: Clinical data 
  doc: Clinical phenotype data
  format: http://edamontology.org/data_2526 #text data
  

outputs:  
 gsea_node_table:
  type: File
  outputSource: gsea_analysis/gsea_node_table
  label: GSEA node table
  doc: Node table for GSEA network result visualisation in Cytoscape
  
 gsea_edge_table:
  type: File
  outputSource: gsea_analysis/gsea_edge_table
  label: GSEA edge table   
  doc: Edge table for GSEA network result visualisation in Cytoscape
  
 net_data_table:
  type: File
  outputSource: gsea_analysis/net_data_table
  label: Network data table
  doc: GSEA network data table
  
 add_net_data_table:
  type: File
  outputSource: gsea_analysis/add_net_data_table
  label:  Additional network data table  
  doc: GSEA network additional visualisation table    

 module_eigengene_net_nodetable:
  type: File
  outputSource: ME_network/module_eigengene_net_nodetable
  label: Module eigengene node data  
  doc: Module eigenegene correlation network node table in Cytoscape
  
 module_eigengene_net_edgetable:
  type: File
  outputSource: ME_network/module_eigengene_net_edgetable
  label: Module eigengene edge data
  doc: Module eigenegene correlation network edge table in Cytoscape 
   
 
steps:
  pre_processing:
    run: pre-processing.cwl
    in: 
      raw_counts: raw_counts
      bgxfile: bgxfile
      
    out:
      [norm_counts]
    label: Data pre-processing 
    doc: Pre-processing input data 

  quality_control:
    run: quality_control.cwl
    in:
       norm_counts: pre_processing/norm_counts
       clinical_data: clinical_data
       
    out:
      [qc_norm_counts, qc_clinical]
    label: Quality control
    doc: Quality control of input data 
      
  DEG:
    run: DEG.cwl
    in:
       qc_norm_counts: quality_control/qc_norm_counts
       qc_clinical: quality_control/qc_clinical
       
    out:
      [deg_data]  
    label: Differential gene expression analysis
    doc: Differential gene expression analysis using limma   
      
  wgcna:
    run: wgcna.cwl
    in:
       qc_norm_counts: quality_control/qc_norm_counts
       qc_clinical: quality_control/qc_clinical
       
    out:
      [module_trait, module_eigengene, module_trait_p]
    label: Weighted gene co-expression network analysis
    doc: Weighted gene co-expression network analysis for module identification
    
  gsea_analysis:
    run: gsea_analysis.cwl
    in:
       deg_data: DEG/deg_data
       
    out:
      [gsea_node_table, gsea_edge_table, net_data_table, add_net_data_table] 
    label: Gene set enrichment analysis and visualisation
    doc: Gene set enrichment analysis and visualisation in Cytoscape  
      
  ME_network:
    run: ME-network.cwl
    in:
       module_eigengene: wgcna/module_eigengene
       module_trait_p: wgcna/module_trait_p
       
    out:
      [module_eigengene_net_nodetable, module_eigengene_net_edgetable]     
    label: Module eigenegene correlation network
    doc: Module eigenegene correlation network and visualisation in Cytoscape      
      
  module_ora:
    run: module_ora.cwl
    in:
       module_trait: wgcna/module_trait
       
    out:
      []   
    label: Module enrichment analysis
    doc: Module enrichment analysis with significantly correlating modules                   
      
s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0003-2878-4049
    s:email: mailto:a.iyer@maastrichtuniversity.nl
    s:name: Aishwarya Iyer

#TODO Add citation
s:codeRepository: https://github.com/aish181095/PSVD-transcriptomics-workflow
#TOOD Add license

$namespaces:
 s: https://schema.org/
 edam: http://edamontology.org/
 

$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - https://edamontology.org/EDAM_1.25.owl
          
      
      
      
