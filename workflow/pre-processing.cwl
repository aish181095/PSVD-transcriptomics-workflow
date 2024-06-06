#!/usr/bin/env cwltool

cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: ashiyer/psvd_workflow:r-4.4-packages

requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: pre_processing.R
      entry:
        $include: ../scripts/pre_processing.R

baseCommand: ["Rscript", "pre_processing.R"]

inputs: 
 raw_counts:
  type: File
  inputBinding:
   position: 1
 
 bgxfile:
  type: File
  inputBinding:
   position: 2

outputs: 
 norm_counts:
  type: File
  outputBinding:
   glob: normalised.data.txt
