#!/usr/bin/env cwltool

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: ashiyer/psvd_workflow:r-4.4-packages
    
requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: ME-network.R
      entry:
        $include: ../scripts/ME-network.R

baseCommand: ["Rscript", "ME-network.R"]

inputs:
 module_eigengene:
   type: File
   inputBinding:
    position: 1
    
 module_trait_p:
   type: File
   inputBinding:
    position: 1      
   
 
outputs: 
 module_eigengene_net_nodetable:
  type: File
  outputBinding:
     glob: corrnet-node-table.txt
      
 module_eigengene_net_edgetable:
  type: File
  outputBinding:
     glob: corrnet-edges-table.txt    
      
  
