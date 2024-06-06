#!/usr/bin/env cwltool

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: ashiyer/psvd_workflow:r-4.4-packages
    
requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: gsea_analysis.R
      entry:
        $include: ../scripts/gsea_analysis.R

baseCommand: ["Rscript", "gsea_analysis.R"]

inputs:
 deg_data:
   type: File
   inputBinding:
    position: 1
   
 
outputs: 
  gsea_node_table:
   type: File
   outputBinding:
      glob: gsea-cytoscape-node-table.txt
      
  gsea_edge_table:
   type: File
   outputBinding:
      glob: gsea-cytoscape-edges-table.txt    
      
  net_data_table:
   type: File
   outputBinding:
      glob: gsea-network-data-table.txt
      
      
  add_net_data_table:
   type: File
   outputBinding:
     glob: add-network-data-table.txt            
  
  
    
  
