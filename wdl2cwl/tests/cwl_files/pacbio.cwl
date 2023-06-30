cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: mergePacBio
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                set -e
                mkdir -p \$(dirname $(inputs.outputPathMergedReport))
                pacbio_merge \
                --reports $(inputs.reports.map(function(el) {return el.path}).join(" ")) \
                --json-output $(inputs.outputPathMergedReport)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/redmar_van_den_berg/pacbio-merge:0.2
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = inputs.memory.match(/[a-zA-Z]+/g).join("");
            var value = parseInt(`${inputs.memory}`.match(/[0-9]+/g));
            var memory = "";
            if(unit==="KiB") memory = value/1024;
            else if(unit==="MiB") memory = value;
            else if(unit==="GiB") memory = value*1024;
            else if(unit==="TiB") memory = value*1024*1024;
            else if(unit==="B") memory = value/(1024*1024);
            else if(unit==="KB" || unit==="K") memory = (value*1000)/(1024*1024);
            else if(unit==="MB" || unit==="M") memory = (value*(1000*1000))/(1024*1024);
            else if(unit==="GB" || unit==="G") memory = (value*(1000*1000*1000))/(1024*1024);
            else if(unit==="TB" || unit==="T") memory = (value*(1000*1000*1000*1000))/(1024*1024);
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
    inputs:
      - id: reports
        doc: The PacBio report files to merge.
        type:
            items: File
            type: array
      - id: outputPathMergedReport
        doc: The location the merged PacBio report file should be written to.
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 4G
        type: string
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/redmar_van_den_berg/pacbio-merge:0.2
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.reports.length == 0) {throw "reports must contain
            at least one item.";} else { return "";}}
    outputs:
      - id: outputMergedReport
        doc: The PacBio reports merged into one.
        type: File
        outputBinding:
            glob: $(inputs.outputPathMergedReport)
  - cwlVersion: v1.2
    id: ccsChunks
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                set -e
                python <<CODE
                for i in range(1, $(inputs.chunkCount) + 1):
                    print(i, $(inputs.chunkCount), sep="/")
                CODE
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: python:3.7-slim
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = inputs.memory.match(/[a-zA-Z]+/g).join("");
            var value = parseInt(`${inputs.memory}`.match(/[0-9]+/g));
            var memory = "";
            if(unit==="KiB") memory = value/1024;
            else if(unit==="MiB") memory = value;
            else if(unit==="GiB") memory = value*1024;
            else if(unit==="TiB") memory = value*1024*1024;
            else if(unit==="B") memory = value/(1024*1024);
            else if(unit==="KB" || unit==="K") memory = (value*1000)/(1024*1024);
            else if(unit==="MB" || unit==="M") memory = (value*(1000*1000))/(1024*1024);
            else if(unit==="GB" || unit==="G") memory = (value*(1000*1000*1000))/(1024*1024);
            else if(unit==="TB" || unit==="T") memory = (value*(1000*1000*1000*1000))/(1024*1024);
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
    inputs:
      - id: chunkCount
        doc: The number of chunks to create.
        type: int
      - id: memory
        doc: The amount of memory this job will use.
        default: 4G
        type: string
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: python:3.7-slim
        type: string
    baseCommand:
      - bash
      - script.bash
    stdout: _stdout
    outputs:
      - id: chunks
        doc: The chunks created based on `chunkCount`.
        type:
            items: string
            type: array
        outputBinding:
            loadContents: true
            glob: _stdout
            outputEval: |-
                ${
                  var contents = self[0].contents.replace(/\r\n$/, "").replace(/\n$/, "").replace(/\r$/, "");
                  // ^ remove any trailing newline to prevent a null being returned
                  return contents.split(/\r\n|\r|\n/);
                }
