cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: Build
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                centrifuge-build \
                --threads $(inputs.threads) \
                $(inputs.disableDifferenceCover ? "--nodc" : "") \
                $(inputs.offrate === null ? "" : "--offrate " + inputs.offrate) \
                $(inputs.ftabChars === null ? "" : "--ftabchars " + inputs.ftabChars) \
                $(inputs.kmerCount === null ? "" : "--kmer-count " + inputs.kmerCount) \
                $(inputs.sizeTable === null ? "" : "--size-table " + inputs.sizeTable.path) \
                --conversion-table $(inputs.conversionTable.path) \
                --taxonomy-tree $(inputs.taxonomyTree.path) \
                --name-table $(inputs.nameTable.path) \
                $(inputs.referenceFile.path) \
                $("/" + inputs.outputPrefix + inputs.indexBasename)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5
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
      - id: disableDifferenceCover
        doc: Disable use of the difference-cover sample.
        default: false
        type: boolean
      - id: conversionTable
        doc: List of UIDs (unique ID) and corresponding taxonomic IDs.
        type: File
      - id: taxonomyTree
        doc: Taxonomic tree (e.g. nodes.dmp).
        type: File
      - id: nameTable
        doc: Name table (e.g. names.dmp).
        type: File
      - id: referenceFile
        doc: A comma-separated list of fasta files containing the reference sequences
            to be aligned to.
        type: File
      - id: indexBasename
        doc: The basename of the index files to write.
        default: centrifuge_index
        type: string
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: offrate
        doc: The number of rows marked by the indexer.
        type:
          - int
          - 'null'
      - id: ftabChars
        doc: Calculate an initial BW range with respect to this character.
        type:
          - int
          - 'null'
      - id: kmerCount
        doc: Use <int> as kmer-size for counting the distinct number of k-mers in
            the input sequences.
        type:
          - int
          - 'null'
      - id: sizeTable
        doc: List of taxonomic IDs and lengths of the sequences belonging to the same
            taxonomic IDs.
        type:
          - File
          - 'null'
      - id: threads
        doc: The number of threads to be used.
        default: 5
        type: int
      - id: memory
        doc: The amount of memory available to the job.
        default: 20G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 2880
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: index
        doc: Generated centrifuge index.
        type:
            items: File
            type: array
        outputBinding:
            glob: $(inputs.outputPrefix + "/" + inputs.indexBasename + "*.cf")
  - cwlVersion: v1.2
    id: Classify
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                indexBasename="\$(basename $(inputs.indexFiles[0].replace("\\.[0-9]\\.cf", "") ))"
                for file in $(inputs.indexFiles.map(function(el) {return el.path}).join(" "))
                do
                    ln ${file} $PWD/"\$(basename ${file})"
                done
                centrifuge \
                $({ "fastq": "-q", "fasta": "-f", "qseq": "--qseq", "raw": "-r", "sequences": "-c" }[inputs.inputFormat]) \
                $(inputs.phred64 ? "--phred64" : "--phred33") \
                --min-hitlen $(inputs.minHitLength) \
                --threads $(inputs.threads) \
                $(inputs.trim5 === null ? "" : "--trim5 " + inputs.trim5) \
                $(inputs.trim3 === null ? "" : "--trim3 " + inputs.trim3) \
                $(inputs.reportMaxDistinct === null ? "" : "-k " + inputs.reportMaxDistinct) \
                $(inputs.hostTaxIDs === null ? "" : "--host-taxids " + inputs.hostTaxIDs) \
                $(inputs.excludeTaxIDs === null ? "" : "--exclude-taxids " + inputs.excludeTaxIDs) \
                -x $PWD/${indexBasename} \
                $(inputs.read2.length > 0 ? "-1" : "-U") $(inputs.read1.map(function(el) {return el.path}).join(",")) \
                $(inputs.read2.length > 0 ? "-2" : "") $(inputs.read2.map(function(el) {return el.path}).join(",")) \
                $("-S " + inputs.outputPrefix + "_classification.tsv") \
                $("--report-file " + inputs.outputPrefix + "_output_report.tsv")
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5
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
      - id: read1
        doc: List of files containing mate 1s, or unpaired reads.
        type:
            items: File
            type: array
      - id: read2
        doc: List of files containing mate 2s.
        default: []
        type:
            items: File
            type: array
      - id: inputFormat
        doc: The format of the read file(s).
        default: fastq
        type: string
      - id: phred64
        doc: If set to true, phred+64 encoding is used.
        default: false
        type: boolean
      - id: minHitLength
        doc: Minimum length of partial hits.
        default: 22
        type: int
      - id: indexFiles
        doc: The files of the index for the reference genomes.
        type:
            items: File
            type: array
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: trim5
        doc: Trim <int> bases from 5' (left) end of each read before alignment.
        type:
          - int
          - 'null'
      - id: trim3
        doc: Trim <int> bases from 3' (right) end of each read before alignment.
        type:
          - int
          - 'null'
      - id: reportMaxDistinct
        doc: It searches for at most <int> distinct, primary assignments for each
            read or pair.
        type:
          - int
          - 'null'
      - id: hostTaxIDs
        doc: A comma-separated list of taxonomic IDs that will be preferred in classification
            procedure.
        type:
          - string
          - 'null'
      - id: excludeTaxIDs
        doc: A comma-separated list of taxonomic IDs that will be excluded in classification
            procedure.
        type:
          - string
          - 'null'
      - id: threads
        doc: The number of threads to be used.
        default: 4
        type: int
      - id: memory
        doc: The amount of memory available to the job.
        default: 16G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 2880
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.read1.length == 0) {throw "read1 must contain at least
            one item.";} else { return "";}}
      - valueFrom: ${if (inputs.indexFiles.length == 0) {throw "indexFiles must contain
            at least one item.";} else { return "";}}
    outputs:
      - id: classification
        doc: File with the classification results.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_classification.tsv")
      - id: report
        doc: File with a classification summary.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_output_report.tsv")
  - cwlVersion: v1.2
    id: Inspect
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                indexBasename="\$(basename $(inputs.indexFiles[0].replace("\\.[0-9]\\.cf", "") ))"
                for file in $(inputs.indexFiles.map(function(el) {return el.path}).join(" "))
                do
                    ln ${file} $PWD/"\$(basename ${file})"
                done
                centrifuge-inspect \
                $({ "fasta": "", "names": "--names", "summary": "--summary", "conversionTable": "--conversion-table", "taxonomyTree": "--taxonomy-tree", "nameTable": "--name-table", "sizeTable": "--size-table" }[inputs.printOption]) \
                $(inputs.across === null ? "" : "--across " + inputs.across) \
                $PWD/${indexBasename} \
                > $("/" + inputs.outputPrefix + inputs.printOption)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5
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
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: printOption
        doc: The output option for inspect (fasta, summary, conversionTable, taxonomyTree,
            nameTable, sizeTable)
        default: fasta
        type: string
      - id: indexFiles
        doc: The files of the index for the reference genomes.
        type:
            items: File
            type: array
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: across
        doc: When printing fasta output, output a newline character every <int> bases.
        type:
          - int
          - 'null'
      - id: memory
        doc: The amount of memory available to the job.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 1
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.indexFiles.length == 0) {throw "indexFiles must contain
            at least one item.";} else { return "";}}
    outputs:
      - id: inspectResult
        doc: Output file according to output option.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "/" + inputs.printOption)
  - cwlVersion: v1.2
    id: KReport
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                indexBasename="\$(basename $(inputs.indexFiles[0].replace("\\.[0-9]\\.cf", "") ))"
                for file in $(inputs.indexFiles.map(function(el) {return el.path}).join(" "))
                do
                    ln ${file} $PWD/"\$(basename ${file})"
                done
                centrifuge-kreport \
                -x $PWD/${indexBasename} \
                $(inputs.noLCA ? "--no-lca" : "") \
                $(inputs.showZeros ? "--show-zeros" : "") \
                $(inputs.isCountTable ? "--is-count-table" : "") \
                $(inputs.minimumScore === null ? "" : "--min-score " + inputs.minimumScore) \
                $(inputs.minimumLength === null ? "" : "--min-length " + inputs.minimumLength) \
                $(inputs.classification.path) \
                > $("_kreport.tsv" + inputs.outputPrefix)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5
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
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: classification
        doc: File with centrifuge classification results.
        type: File
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: indexFiles
        doc: The files of the index for the reference genomes.
        type:
            items: File
            type: array
      - id: noLCA
        doc: Do not report the lca of multiple assignments, but report count fractions
            at the taxa.
        default: false
        type: boolean
      - id: showZeros
        doc: Show clades that have zero reads.
        default: false
        type: boolean
      - id: isCountTable
        doc: The format of the file is taxID<tab>COUNT.
        default: false
        type: boolean
      - id: minimumScore
        doc: Require a minimum score for reads to be counted.
        type:
          - int
          - 'null'
      - id: minimumLength
        doc: Require a minimum alignment length to the read.
        type:
          - int
          - 'null'
      - id: memory
        doc: The amount of memory available to the job.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 10
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.indexFiles.length == 0) {throw "indexFiles must contain
            at least one item.";} else { return "";}}
    outputs:
      - id: KrakenReport
        doc: File with kraken style report.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_kreport.tsv")
  - cwlVersion: v1.2
    id: KTimportTaxonomy
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                cat $(inputs.inputFile.path) | cut -f 1,3 > kronaInput.krona
                ktImportTaxonomy kronaInput.krona
                cp taxonomy.krona.html $("_krona.html" + inputs.outputPrefix)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/krona:v2.7.1_cv1
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
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: inputFile
        doc: File with centrifuge classification results.
        type: File
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: memory
        doc: The amount of memory available to the job.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 1
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/krona:v2.7.1_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: kronaPlot
        doc: Krona taxonomy plot html file.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_krona.html")
