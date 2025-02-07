cwlVersion: v1.1
class: CommandLineTool
label: OME-TIFF metadata normalization via bftools
requirements:
  DockerRequirement:
    dockerPull: hubmap/epic-segmentation-mask-convert:latest
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
  calculated_metadata_json:
    type: File
    outputBinding:
      glob: calculated_metadata.json
