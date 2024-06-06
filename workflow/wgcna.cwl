#!/usr/bin/env cwltool

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: ashiyer/psvd_workflow:r-4.4-packages

requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: wgcna.R
      entry:
        $include: ../scripts/wgcna.R

baseCommand: ["Rscript", "wgcna.R"]

inputs:
 qc_norm_counts:
   type: File
   inputBinding:
    position: 1
   
 qc_clinical:
   type: File
   inputBinding:
    position: 2   
   

outputs: 
  module_trait:
   type: File
   outputBinding:
      glob: module-trait.txt
     
  module_eigengene:
   type: File
   outputBinding:
      glob: Module-eigengene.txt 
      
  module_trait_p:
   type: File
   outputBinding:
      glob: moduleTraitPvalue.txt
