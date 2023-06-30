cwlVersion: v1.2
id: ValidateLoom
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |2

            cacheInvalidationRandomString=4

            echo Starting checksum generation...
            calculated_loom_file_checksum=\$( md5sum < $(inputs.loom_file === null ? "" : inputs.loom_file.path) | awk '{print $1}' )
            echo Checksum generation complete

            if [ "$calculated_loom_file_checksum" == "$(inputs.expected_loom_file_checksum)" ]
            then
                echo Computed and expected loom file hashes match \( "$calculated_loom_file_checksum" \)
            printf PASS > result.txt
            else
                echo Computed \( $calculated_loom_file_checksum \) and expected \( $(inputs.expected_loom_file_checksum) \) loom file hashes do not match
               printf FAIL > result.txt
            fi
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: ubuntu:16.04
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 3576.2786865234375
    outdirMin: '$((Math.ceil((function(size_of=0){inputs.loom_file === null ? "" :
        inputs.loom_file.path.forEach(function(element){ if (element) {size_of +=
        element.size}})}) / 1000^3 * 1.1) ) * 1024)'
inputs:
  - id: loom_file
    type:
      - File
      - 'null'
  - id: expected_loom_file_checksum
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: result
    type: string
    outputBinding:
        loadContents: true
        glob: result.txt
        outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
