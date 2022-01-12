class: CommandLineTool
id: ValidateLoom
inputs:
  - id: loom_file
    type:
      - File
      - 'null'
  - id: expected_loom_file_checksum
    type: string
outputs:
  - id: result
    type: string
    outputBinding:
        loadContents: true
        glob: result.txt
        outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
requirements:
  - class: DockerRequirement
    dockerPull: ubuntu:16.04
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

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
  - class: ResourceRequirement
    coresMin: 1
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
