cwlVersion: v1.2
id: Lima
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |2

            set -e
            mkdir -p "\$(dirname $(inputs.outputPrefix))"
            lima \
            $({ "same": "--same", "different": "--different", "neighbors": "--neighbors" }[inputs.libraryDesign]) \
            $(inputs.scoreFullPass ? "--score-full-pass" : "") \
            --max-scored-barcode-pairs $(inputs.maxScoredBarcodePairs) \
            --max-scored-barcodes $(inputs.maxScoredBarcodes) \
            --max-scored-adapters $(inputs.maxScoredAdapters) \
            --min-passes $(inputs.minPasses) \
            --min-length $(inputs.minLength) \
            --max-input-length $(inputs.maxInputLength) \
            --min-ref-span $(inputs.minRefSpan) \
            --min-scoring-regions $(inputs.minScoringRegion) \
            --min-score $(inputs.minScore) \
            --min-end-score $(inputs.minEndScore) \
            --min-signal-increase $(inputs.minSignalIncrease) \
            --min-score-lead $(inputs.minScoreLead) \
            $(inputs.ccsMode ? "--ccs" : "") \
            $(inputs.splitBamNamed ? "--split-bam-named" : "") \
            --scored-adapter-ratio $(inputs.scoredAdapterRatio) \
            --peek $(inputs.peek) \
            --guess $(inputs.guess) \
            --guess-min-count $(inputs.guessMinCount) \
            $(inputs.peekGuess ? "--peek-guess" : "") \
            --log-level $(inputs.logLevel) \
            --num-threads $(inputs.threads) \
            $("--log-file " + inputs.outputPrefix + ".lima.stderr.log") \
            $(inputs.inputBamFile.path) \
            $(inputs.barcodeFile.path) \
            $(".bam" + inputs.outputPrefix)

            dirName="\$(dirname $(inputs.outputPrefix))"
            find "\$(cd ${dirName}; pwd)" -name "*.bam" > bamFiles.txt
            find "\$(cd ${dirName}; pwd)" -name "*.bam.pbi" > bamIndexes.txt
            find "\$(cd ${dirName}; pwd)" -name "*.consensusreadset.xml" > consensusreadset.txt
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/lima:2.2.0--h9ee0642_0
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
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
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes * 60)
inputs:
  - id: libraryDesign
    doc: Barcode structure of the library design.
    default: same
    type: string
  - id: scoreFullPass
    doc: Only use subreads flanked by adapters for barcode identification.
    default: false
    type: boolean
  - id: maxScoredBarcodePairs
    doc: Only use up to N barcode pair regions to find the barcode, 0 means use all.
    default: 0
    type: int
  - id: maxScoredBarcodes
    doc: Analyze at maximum the provided number of barcodes per ZMW, 0 means deactivated.
    default: 0
    type: int
  - id: maxScoredAdapters
    doc: Analyze at maximum the provided number of adapters per ZMW, 0 means deactivated.
    default: 0
    type: int
  - id: minPasses
    doc: Minimal number of full passes.
    default: 0
    type: int
  - id: minLength
    doc: Minimum sequence length after clipping.
    default: 50
    type: int
  - id: maxInputLength
    doc: Maximum input sequence length, 0 means deactivated.
    default: 0
    type: int
  - id: minRefSpan
    doc: Minimum reference span relative to the barcode length.
    default: 0.5
    type: float
  - id: minScoringRegion
    doc: Minimum number of barcode regions with sufficient relative span to the barcode
        length.
    default: 1
    type: int
  - id: minScore
    doc: Reads below the minimum barcode score are removed from downstream analysis.
    default: 0
    type: int
  - id: minEndScore
    doc: Minimum end barcode score threshold is applied to the individual leading
        and trailing ends.
    default: 0
    type: int
  - id: minSignalIncrease
    doc: The minimal score difference, between first and combined, required to call
        a barcode pair different.
    default: 10
    type: int
  - id: minScoreLead
    doc: The minimal score lead required to call a barcode pair significant.
    default: 10
    type: int
  - id: ccsMode
    doc: Ccs mode, use optimal alignment options.
    default: false
    type: boolean
  - id: splitBamNamed
    doc: Split bam output by resolved barcode pair name.
    default: false
    type: boolean
  - id: scoredAdapterRatio
    doc: Minimum ratio of scored vs sequenced adapters.
    default: 0.25
    type: float
  - id: peek
    doc: Demux the first N ZMWs and return the mean score, 0 means peeking deactivated.
    default: 0
    type: int
  - id: guess
    doc: Try to guess the used barcodes, using the provided mean score threshold,
        0 means guessing deactivated.
    default: 0
    type: int
  - id: guessMinCount
    doc: Minimum number of ZMWs observed to whitelist barcodes.
    default: 0
    type: int
  - id: peekGuess
    doc: Try to infer the used barcodes subset, by peeking at the first 50,000 ZMWs.
    default: false
    type: boolean
  - id: logLevel
    doc: 'Set log level. Valid choices: (TRACE, DEBUG, INFO, WARN, FATAL).'
    default: WARN
    type: string
  - id: inputBamFile
    doc: Bam input file.
    type: File
  - id: barcodeFile
    doc: Barcode/primer fasta file.
    type: File
  - id: outputPrefix
    doc: Output directory path + output file prefix.
    type: string
  - id: threads
    doc: The number of threads to be used.
    default: 2
    type: int
  - id: memory
    doc: The amount of memory available to the job.
    default: 2G
    type: string
  - id: timeMinutes
    doc: The maximum amount of time the job will run in minutes.
    default: 30
    type: int
  - id: dockerImage
    doc: The docker image used for this task. Changing this may result in errors which
        the developers may choose not to address.
    default: quay.io/biocontainers/lima:2.2.0--h9ee0642_0
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: limaBam
    doc: Demultiplexed reads output file(s).
    type:
        items: File
        type: array
    outputBinding:
        loadContents: true
        glob: bamFiles.txt
        outputEval: |-
            ${
              var contents = self[0].contents.replace(/\r\n$/, "").replace(/\n$/, "").replace(/\r$/, "");
              // ^ remove any trailing newline to prevent a null being returned
              return contents.split(/\r\n|\r|\n/);
            }
  - id: limaBamIndex
    doc: Index of demultiplexed reads output file(s).
    type:
        items: File
        type: array
    outputBinding:
        loadContents: true
        glob: bamIndexes.txt
        outputEval: |-
            ${
              var contents = self[0].contents.replace(/\r\n$/, "").replace(/\n$/, "").replace(/\r$/, "");
              // ^ remove any trailing newline to prevent a null being returned
              return contents.split(/\r\n|\r|\n/);
            }
  - id: limaXml
    doc: Xml file of the subreadset(s).
    type:
        items: File
        type: array
    outputBinding:
        loadContents: true
        glob: consensusreadset.txt
        outputEval: |-
            ${
              var contents = self[0].contents.replace(/\r\n$/, "").replace(/\n$/, "").replace(/\r$/, "");
              // ^ remove any trailing newline to prevent a null being returned
              return contents.split(/\r\n|\r|\n/);
            }
  - id: limaStderr
    doc: Lima stderr log file.
    type: File
    outputBinding:
        glob: $(inputs.outputPrefix + ".lima.stderr.log")
  - id: limaJson
    doc: Lima json file.
    type: File
    outputBinding:
        glob: $(inputs.outputPrefix + ".json")
  - id: limaCounts
    doc: Lima counts file.
    type: File
    outputBinding:
        glob: $(inputs.outputPrefix + ".lima.counts")
  - id: limaReport
    doc: Lima report file.
    type: File
    outputBinding:
        glob: $(inputs.outputPrefix + ".lima.report")
  - id: limaSummary
    doc: Lima summary file.
    type: File
    outputBinding:
        glob: $(inputs.outputPrefix + ".lima.summary")
