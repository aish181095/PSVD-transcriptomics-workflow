#!/usr/bin/env cwltool

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: ashiyer/psvd_workflow:r-4.4-packages

requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: quality_control.R
      entry:
        $include: ../scripts/quality_control.R

baseCommand: ["Rscript", "quality_control.R"]

inputs:
 norm_counts:
   type: File
   inputBinding:
    position: 1
   
 clinical_data:
  type: File
  inputBinding:
   position: 2   
   

outputs: 
  qc_norm_counts:
   type: File
   outputBinding:
      glob: normalised-data-afteroutlier.txt
     
  qc_clinical:
   type: File
   outputBinding:
      glob: clinical-data-afteroutlier.txt    
      
  
    
  
