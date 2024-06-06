#!/usr/bin/env cwltool

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: ashiyer/psvd_workflow:r-4.4-packages

requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: module-ora.R
      entry:
        $include: ../scripts/module-ora.R

baseCommand: ["Rscript", "module-ora.R"]
inputs:
 module_trait:
   type: File
   inputBinding:
    position: 1
   
    

outputs: 
  []
