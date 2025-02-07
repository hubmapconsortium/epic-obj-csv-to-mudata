cwlVersion: v1.1
class: CommandLineTool
label: OME-TIFF metadata normalization via bftools
requirements:
  DockerRequirement:
    dockerPull: hubmap/epic-segmentation-mask-convert:0.3.1
  InlineJavascriptRequirement: {}
baseCommand: /opt/obj_feature_to_mudata.py

inputs:
  data_dir:
    type: Directory
    inputBinding:
      position: 0
outputs:
  mudata_file:
    type: File
    outputBinding:
      glob: objects.h5mu
  metadata_json:
    type: File
    outputBinding:
      glob: metadata.json
