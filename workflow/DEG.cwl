#!/usr/bin/env cwltool

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: ashiyer/psvd_workflow:r-4.4-packages

requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: DEG.R
      entry:
        $include: ../scripts/DEG.R
        
baseCommand: ["Rscript", "DEG.R"]

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
 deg_data:
  type: File
  outputBinding:
   glob: deg.data.txt        
