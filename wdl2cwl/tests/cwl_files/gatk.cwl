cwlVersion: v1.2
$graph:
  - class: CommandLineTool
    id: AnnotateIntervals
    inputs:
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: annotatedIntervalsPath
        default: intervals.annotated.tsv
        type: string
      - id: intervals
        type: File
      - id: intervalMergingRule
        default: OVERLAPPING_ONLY
        type: string
      - id: featureQueryLookahead
        default: 1000000
        type: int
      - id: mappabilityTrack
        type:
          - File
          - 'null'
      - id: segmentalDuplicationTrack
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 2G
        type: string
      - id: memory
        default: 3G
        type: string
      - id: timeMinutes
        default: 5
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: annotatedIntervals
        type: File
        outputBinding:
            glob: $(inputs.annotatedIntervalsPath)
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.annotatedIntervalsPath))"
                mkdir reference_dir
                ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
                ln -s $(inputs.referenceFastaDict.path) reference_dir/$(inputs.referenceFastaDict.basename)
                ln -s $(inputs.referenceFastaFai.path) reference_dir/$(inputs.referenceFastaFai.basename)
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                AnnotateIntervals \
                -R reference_dir/$(inputs.referenceFasta.basename) \
                -L $(inputs.intervals.path) \
                $(inputs.mappabilityTrack === null ? "" : "--mappability-track  " + inputs.mappabilityTrack.path) \
                $(inputs.segmentalDuplicationTrack === null ? "" : "--segmental-duplication-track " + inputs.segmentalDuplicationTrack.path) \
                --feature-query-lookahead $(inputs.featureQueryLookahead) \
                --interval-merging-rule $(inputs.intervalMergingRule) \
                -O $(inputs.annotatedIntervalsPath)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: ApplyBQSR
    inputs:
      - id: inputBam
        type: File
      - id: inputBamIndex
        type: File
      - id: outputBamPath
        type: string
      - id: recalibrationReport
        type: File
      - id: sequenceGroupInterval
        default: []
        type:
            items: File
            type: array
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: javaXmxMb
        default: 2048
        type: int
      - id: memoryMb
        type:
          - int
          - 'null'
      - id: timeMinutes
        default: 120
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: recalibratedBam
        type: File
        outputBinding:
            glob: $(inputs.outputBamPath)
      - id: recalibratedBamIndex
        type: File
        outputBinding:
            glob: $(inputs.outputBamPath.replace("\.bam$", ".bai") )
      - id: recalibratedBamMd5
        type: File
        outputBinding:
            glob: $(inputs.outputBamPath + ".md5")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputBamPath))"
                mkdir bam_dir
                ln -s $(inputs.inputBam.path) bam_dir/$(inputs.inputBam.basename)
                ln -s $(inputs.inputBamIndex.path) bam_dir/$(inputs.inputBamIndex.basename)
                mkdir reference_dir
                ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
                ln -s $(inputs.referenceFastaDict.path) reference_dir/$(inputs.referenceFastaDict.basename)
                ln -s $(inputs.referenceFastaFai.path) reference_dir/$(inputs.referenceFastaFai.basename)
                gatk --java-options '-Xmx$(inputs.javaXmxMb)M -XX:ParallelGCThreads=1' \
                ApplyBQSR \
                --create-output-bam-md5 \
                --add-output-sam-program-record \
                -R reference_dir/$(inputs.referenceFasta.basename) \
                -I bam_dir/$(inputs.inputBam.basename) \
                --use-original-qualities \
                -O $(inputs.outputBamPath) \
                -bqsr $(inputs.recalibrationReport.path) \
                --static-quantized-quals 10 \
                --static-quantized-quals 20 \
                --static-quantized-quals 30 \
                $(inputs.sequenceGroupInterval.length > 0 ? "-L" : "") $(inputs.sequenceGroupInterval.map(function(el) {return el.path}).join(" -L "))
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = "M";
            var value = parseInt(`${inputs.javaXmxMb + 512}`.match(/[0-9]+/g));
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: BaseRecalibrator
    inputs:
      - id: inputBam
        type: File
      - id: inputBamIndex
        type: File
      - id: recalibrationReportPath
        type: string
      - id: sequenceGroupInterval
        default: []
        type:
            items: File
            type: array
      - id: knownIndelsSitesVCFs
        default: []
        type:
            items: File
            type: array
      - id: knownIndelsSitesVCFIndexes
        default: []
        type:
            items: File
            type: array
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: dbsnpVCF
        type:
          - File
          - 'null'
      - id: dbsnpVCFIndex
        type:
          - File
          - 'null'
      - id: javaXmxMb
        default: 1024
        type: int
      - id: memoryMb
        type:
          - int
          - 'null'
      - id: timeMinutes
        default: 120
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: recalibrationReport
        type: File
        outputBinding:
            glob: $(inputs.recalibrationReportPath)
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.recalibrationReportPath))"
                mkdir bam_dir
                ln -s $(inputs.inputBam.path) bam_dir/$(inputs.inputBam.basename)
                ln -s $(inputs.inputBamIndex.path) bam_dir/$(inputs.inputBamIndex.basename)
                mkdir reference_dir
                ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
                ln -s $(inputs.referenceFastaDict.path) reference_dir/$(inputs.referenceFastaDict.basename)
                ln -s $(inputs.referenceFastaFai.path) reference_dir/$(inputs.referenceFastaFai.basename)
                gatk --java-options '-Xmx$(inputs.javaXmxMb)M -XX:ParallelGCThreads=1' \
                BaseRecalibrator \
                -R reference_dir/$(inputs.referenceFasta.basename) \
                -I bam_dir/$(inputs.inputBam.basename) \
                --use-original-qualities \
                -O $(inputs.recalibrationReportPath) \
                $(inputs.knownIndelsSitesVCFs.length > 0 ? "--known-sites" : "") $(inputs.knownIndelsSitesVCFs.map(function(el) {return el.path}).join(" --known-sites ")) \
                $(inputs.dbsnpVCF === null ? "" : "--known-sites " + inputs.dbsnpVCF.path) \
                $(inputs.sequenceGroupInterval.length > 0 ? "-L" : "") $(inputs.sequenceGroupInterval.map(function(el) {return el.path}).join(" -L "))
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = "M";
            var value = parseInt(`${inputs.javaXmxMb + 512}`.match(/[0-9]+/g));
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: CalculateContamination
    inputs:
      - id: tumorPileups
        type: File
      - id: normalPileups
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 12G
        type: string
      - id: memory
        default: 13G
        type: string
      - id: timeMinutes
        default: 180
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: contaminationTable
        type: File
        outputBinding:
            glob: contamination.table
      - id: mafTumorSegments
        type: File
        outputBinding:
            glob: segments.table
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                CalculateContamination \
                -I $(inputs.tumorPileups.path) \
                $(inputs.normalPileups === null ? "" : "-matched " + inputs.normalPileups.path) \
                -O "contamination.table" \
                --tumor-segmentation "segments.table"
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: CallCopyRatioSegments
    inputs:
      - id: outputPrefix
        type: string
      - id: copyRatioSegments
        type: File
      - id: javaXmx
        default: 2G
        type: string
      - id: memory
        default: 3G
        type: string
      - id: timeMinutes
        default: 2
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: calledSegments
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".called.seg")
      - id: calledSegmentsIgv
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".called.igv.seg")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                CallCopyRatioSegments \
                -I $(inputs.copyRatioSegments.path) \
                -O $(inputs.outputPrefix).called.seg
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: CollectAllelicCounts
    inputs:
      - id: allelicCountsPath
        default: allelic_counts.tsv
        type: string
      - id: commonVariantSites
        type: File
      - id: inputBam
        type: File
      - id: inputBamIndex
        type: File
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: commonVariantSitesIndex
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 10G
        type: string
      - id: memory
        default: 11G
        type: string
      - id: timeMinutes
        default: 120
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: allelicCounts
        type: File
        outputBinding:
            glob: $(inputs.allelicCountsPath)
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.allelicCountsPath))"
                mkdir bam_dir
                ln -s $(inputs.inputBam.path) bam_dir/$(inputs.inputBam.basename)
                ln -s $(inputs.inputBamIndex.path) bam_dir/$(inputs.inputBamIndex.basename)
                mkdir reference_dir
                ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
                ln -s $(inputs.referenceFastaDict.path) reference_dir/$(inputs.referenceFastaDict.basename)
                ln -s $(inputs.referenceFastaFai.path) reference_dir/$(inputs.referenceFastaFai.basename)
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                CollectAllelicCounts \
                -L $(inputs.commonVariantSites.path) \
                -I bam_dir/$(inputs.inputBam.basename) \
                -R reference_dir/$(inputs.referenceFasta.basename) \
                -O $(inputs.allelicCountsPath)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: CollectReadCounts
    inputs:
      - id: countsPath
        default: readcounts.hdf5
        type: string
      - id: intervals
        type: File
      - id: inputBam
        type: File
      - id: inputBamIndex
        type: File
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: intervalMergingRule
        default: OVERLAPPING_ONLY
        type: string
      - id: javaXmx
        default: 7G
        type: string
      - id: memory
        default: 8G
        type: string
      - id: timeMinutes
        type:
          - int
          - 'null'
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: counts
        type: File
        outputBinding:
            glob: $(inputs.countsPath)
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.countsPath))"
                mkdir bam_dir
                ln -s $(inputs.inputBam.path) bam_dir/$(inputs.inputBam.basename)
                ln -s $(inputs.inputBamIndex.path) bam_dir/$(inputs.inputBamIndex.basename)
                mkdir reference_dir
                ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
                ln -s $(inputs.referenceFastaDict.path) reference_dir/$(inputs.referenceFastaDict.basename)
                ln -s $(inputs.referenceFastaFai.path) reference_dir/$(inputs.referenceFastaFai.basename)
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                CollectReadCounts \
                -L $(inputs.intervals.path) \
                -I bam_dir/$(inputs.inputBam.basename) \
                -R reference_dir/$(inputs.referenceFasta.basename) \
                --format HDF5 \
                --interval-merging-rule $(inputs.intervalMergingRule) \
                -O $(inputs.countsPath)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.inputBam.path.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3 * 5)  * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: CombineGVCFs
    inputs:
      - id: gvcfFiles
        type:
            items: File
            type: array
      - id: gvcfFilesIndex
        type:
            items: File
            type: array
      - id: intervals
        default: []
        type:
            items: File
            type: array
      - id: outputPath
        type: string
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: javaXmx
        default: 4G
        type: string
      - id: memory
        default: 5G
        type: string
      - id: timeMinutes
        type:
          - int
          - 'null'
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: outputVcf
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVcfIndex
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPath))"
                mkdir wd
                for FILE in $(inputs.gvcfFiles.map(function(el) {return el.path}).join(" ")); do ln -s $FILE wd/\$(basename $FILE) ; done
                for FILE in $(inputs.gvcfFilesIndex.map(function(el) {return el.path}).join(" ")); do ln -s $FILE wd/\$(basename $FILE) ; done
                mkdir reference_dir
                ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
                ln -s $(inputs.referenceFastaDict.path) reference_dir/$(inputs.referenceFastaDict.basename)
                ln -s $(inputs.referenceFastaFai.path) reference_dir/$(inputs.referenceFastaFai.basename)
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                CombineGVCFs \
                -R reference_dir/$(inputs.referenceFasta.basename) \
                -O $(inputs.outputPath) \
                (for FILE in $(inputs.gvcfFiles.map(function(el) {return el.path}).join(" ")); do echo -- "-V wd/"\$(basename $FILE); done) \
                $(inputs.intervals.length > 0 ? "-L" : "") $(inputs.intervals.map(function(el) {return el.path}).join(" -L "))
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.gvcfFiles.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3 * 8)  * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: CombineVariants
    inputs:
      - id: referenceFasta
        type: File
      - id: referenceFastaFai
        type: File
      - id: referenceFastaDict
        type: File
      - id: genotypeMergeOption
        default: UNIQUIFY
        type: string
      - id: filteredRecordsMergeType
        default: KEEP_IF_ANY_UNFILTERED
        type: string
      - id: identifiers
        type:
            items: string
            type: array
      - id: variantVcfs
        type:
            items: File
            type: array
      - id: variantIndexes
        type:
            items: File
            type: array
      - id: outputPath
        type: string
      - id: javaXmx
        default: 12G
        type: string
      - id: memory
        default: 13G
        type: string
      - id: timeMinutes
        default: 180
        type: int
      - id: dockerImage
        default: broadinstitute/gatk3:3.8-1
        type: string
    outputs:
      - id: combinedVcf
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: combinedVcfIndex
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir wd
                for FILE in ${sep(" ", variantVcfs)}; do ln -s $FILE wd/\$(basename $FILE) ; done
                for FILE in ${sep(" ", variantIndexes)}; do ln -s $FILE wd/\$(basename $FILE) ; done
                mkdir reference_dir
                ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
                ln -s $(inputs.referenceFastaDict.path) reference_dir/$(inputs.referenceFastaDict.basename)
                ln -s $(inputs.referenceFastaFai.path) reference_dir/$(inputs.referenceFastaFai.basename)
                mkdir -p "\$(dirname $(inputs.outputPath))"
                # Build "-V:<ID> <file.vcf>" arguments according to IDs
                # and VCFs to merge.
                # Make sure commands are run in bash.
                V_args=\$(bash -c '
                set -eu
                ids=($(inputs.identifiers.join(" ")))
                vars=(\$(for file in ${sep(" ", variantVcfs); do echo wd/\$(basename $file) ; done))
                for (( i = 0; i < ${#ids[@]}; ++i ))
                do
                    printf -- "-V:%s %s " "${ids[i]}" "${vars[i]}"
                done
                ')
                java -Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1 -jar /usr/GenomeAnalysisTK.jar \
                -T CombineVariants \
                -R reference_dir/$(inputs.referenceFasta.basename) \
                --genotypemergeoption $(inputs.genotypeMergeOption) \
                --filteredrecordsmergetype $(inputs.filteredRecordsMergeType) \
                --out $(inputs.outputPath) \
                $V_args
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: broadinstitute/gatk3:3.8-1
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(180 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: CreateReadCountPanelOfNormals
    inputs:
      - id: PONpath
        default: PON.hdf5
        type: string
      - id: readCountsFiles
        type:
            items: File
            type: array
      - id: annotatedIntervals
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 7G
        type: string
      - id: memory
        default: 8G
        type: string
      - id: timeMinutes
        default: 5
        type: int
      - id: dockerImage
        default: broadinstitute/gatk:4.1.8.0
        type: string
    outputs:
      - id: PON
        type: File
        outputBinding:
            glob: $(inputs.PONpath)
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.PONpath))"
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                CreateReadCountPanelOfNormals \
                -I $(inputs.readCountsFiles.map(function(el) {return el.path}).join(" -I ")) \
                $(inputs.annotatedIntervals === null ? "" : "--annotated-intervals " + inputs.annotatedIntervals.path) \
                -O $(inputs.PONpath)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: broadinstitute/gatk:4.1.8.0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(5 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: DenoiseReadCounts
    inputs:
      - id: readCounts
        type: File
      - id: outputPrefix
        type: string
      - id: PON
        type:
          - File
          - 'null'
      - id: annotatedIntervals
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 4G
        type: string
      - id: memory
        default: 5G
        type: string
      - id: timeMinutes
        default: 5
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: standardizedCopyRatios
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".standardizedCR.tsv")
      - id: denoisedCopyRatios
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".denoisedCR.tsv")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                DenoiseReadCounts \
                -I $(inputs.readCounts.path) \
                $(inputs.PON === null ? "" : "--count-panel-of-normals " + inputs.PON.path) \
                $(inputs.annotatedIntervals === null ? "" : "--annotated-intervals " + inputs.annotatedIntervals.path) \
                --standardized-copy-ratios $(inputs.outputPrefix).standardizedCR.tsv \
                --denoised-copy-ratios $(inputs.outputPrefix).denoisedCR.tsv
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(5 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: FilterMutectCalls
    inputs:
      - id: referenceFasta
        type: File
      - id: referenceFastaFai
        type: File
      - id: referenceFastaDict
        type: File
      - id: unfilteredVcf
        type: File
      - id: unfilteredVcfIndex
        type: File
      - id: outputVcf
        type: string
      - id: uniqueAltReadCount
        default: 4
        type: int
      - id: mutect2Stats
        type: File
      - id: contaminationTable
        type:
          - File
          - 'null'
      - id: mafTumorSegments
        type:
          - File
          - 'null'
      - id: artifactPriors
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 12G
        type: string
      - id: memory
        default: 13G
        type: string
      - id: timeMinutes
        default: 60
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: filteredVcf
        type: File
        outputBinding:
            glob: $(inputs.outputVcf)
      - id: filteredVcfIndex
        type: File
        outputBinding:
            glob: $(inputs.outputVcf + ".tbi")
      - id: filteringStats
        type: File
        outputBinding:
            glob: filtering.stats
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputVcf))"
                mkdird unfiltered_vcf_dir
                ln -s $(inputs.unfilteredVcf.path) unfiltered_vcf_dir/$(inputs.unfilteredVcf.basename)
                ln -s $(inputs.unfilteredVcfIndex.path) unfiltered_vcf_dir/$(inputs.unfilteredVcfIndex.basename)
                mkdir reference_dir
                ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
                ln -s $(inputs.referenceFastaDict.path) reference_dir/$(inputs.referenceFastaDict.basename)
                ln -s $(inputs.referenceFastaFai.path) reference_dir/$(inputs.referenceFastaFai.basename)
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                FilterMutectCalls \
                -R reference_dir/$(inputs.referenceFasta.basename) \
                -V unfiltered_vcf_dir/$(inputs.unfilteredVcf.basename) \
                -O $(inputs.outputVcf) \
                $(inputs.contaminationTable === null ? "" : "--contamination-table " + inputs.contaminationTable.path) \
                $(inputs.mafTumorSegments === null ? "" : "--tumor-segmentation " + inputs.mafTumorSegments.path) \
                $(inputs.artifactPriors === null ? "" : "--ob-priors " + inputs.artifactPriors.path) \
                --unique-alt-read-count  $(inputs.uniqueAltReadCount) \
                -stats  $(inputs.mutect2Stats.path) \
                --filtering-stats "filtering.stats" \
                --showHidden
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(60 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: GatherBqsrReports
    inputs:
      - id: inputBQSRreports
        type:
            items: File
            type: array
      - id: outputReportPath
        type: string
      - id: javaXmxMb
        default: 256
        type: int
      - id: memoryMb
        type:
          - int
          - 'null'
      - id: timeMinutes
        default: 1
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: outputBQSRreport
        type: File
        outputBinding:
            glob: $(inputs.outputReportPath)
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputReportPath))"
                gatk --java-options '-Xmx$(inputs.javaXmxMb)M -XX:ParallelGCThreads=1' \
                GatherBQSRReports \
                -I $(inputs.inputBQSRreports.map(function(el) {return el.path}).join(" -I ")) \
                -O $(inputs.outputReportPath)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = "M";
            var value = parseInt(`${256 + inputs.javaXmxMb}`.match(/[0-9]+/g));
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(1 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: GenomicsDBImport
    inputs:
      - id: gvcfFiles
        type:
            items: File
            type: array
      - id: gvcfFilesIndex
        type:
            items: File
            type: array
      - id: intervals
        type:
            items: File
            type: array
      - id: genomicsDBWorkspacePath
        default: genomics_db
        type: string
      - id: genomicsDBTarFile
        default: genomics_db.tar.gz
        type: string
      - id: tmpDir
        type:
          - string
          - 'null'
      - id: javaXmx
        default: 4G
        type: string
      - id: memory
        default: 5G
        type: string
      - id: timeMinutes
        default: 180
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: genomicsDbTarArchive
        type: File
        outputBinding:
            glob: $(inputs.genomicsDBTarFile)
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.genomicsDBWorkspacePath))"
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                GenomicsDBImport \
                -V $(inputs.gvcfFiles.map(function(el) {return el.path}).join(" -V ")) \
                --genomicsdb-workspace-path $(inputs.genomicsDBWorkspacePath) \
                $(inputs.tmpDir === null ? "" : "--tmp-dir " + inputs.tmpDir) \
                -L $(inputs.intervals.map(function(el) {return el.path}).join(" -L "))
                bash -c 'tar -cvzf $(inputs.genomicsDBTarFile) $(inputs.genomicsDBWorkspacePath)/*'
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(180 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: GenotypeGVCFs
    inputs:
      - id: gvcfFile
        type: File
      - id: gvcfFileIndex
        type: File
      - id: outputPath
        type: string
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: annotationGroups
        default:
          - StandardAnnotation
        type:
            items: string
            type: array
      - id: intervals
        type:
          - items: File
            type: array
          - 'null'
      - id: dbsnpVCF
        type:
          - File
          - 'null'
      - id: dbsnpVCFIndex
        type:
          - File
          - 'null'
      - id: pedigree
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 6G
        type: string
      - id: memory
        default: 7G
        type: string
      - id: timeMinutes
        default: 120
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: outputVCF
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVCFIndex
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPath))"
                mkdir wd
                ln -s $(inputs.gvcfFile.path) wd/$(inputs.gvcfFile.basename)
                ln -s $(inputs.gvcfFileIndex.path) wd/$(inputs.gvcfFileIndex.basename)
                mkdir reference_dir
                ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
                ln -s $(inputs.referenceFastaDict.path) reference_dir/$(inputs.referenceFastaDict.basename)
                ln -s $(inputs.referenceFastaFai.path) reference_dir/$(inputs.referenceFastaFai.basename)
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                GenotypeGVCFs \
                -R reference_dir/$(inputs.referenceFasta.basename) \
                -O $(inputs.outputPath) \
                $(inputs.dbsnpVCF === null ? "" : "-D " + inputs.dbsnpVCF.path) \
                $(inputs.pedigree === null ? "" : "--pedigree " + inputs.pedigree.path) \
                $(inputs.annotationGroups.length > 0 ? "-G" : "") $(inputs.annotationGroups.join(" -G ")) \
                -V wd/$(inputs.gvcfFile.basename) \
                $(inputs.intervals === null ? "" : "--only-output-calls-starting-in-intervals") \
                $(inputs.intervals === null ? "" : "-L") $(inputs.intervals.map(function(el) {return el.path}).join(" -L "))
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(120 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: GetPileupSummaries
    inputs:
      - id: sampleBam
        type: File
      - id: sampleBamIndex
        type: File
      - id: variantsForContamination
        type: File
      - id: variantsForContaminationIndex
        type: File
      - id: sitesForContamination
        type: File
      - id: sitesForContaminationIndex
        type: File
      - id: outputPrefix
        type: string
      - id: javaXmx
        default: 12G
        type: string
      - id: memory
        default: 13G
        type: string
      - id: timeMinutes
        default: 120
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: pileups
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "-pileups.table")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir bam_dir
                ln -s $(inputs.sampleBam.path) bam_dir/$(inputs.sampleBam.basename)
                ln -s $(inputs.sampleBamIndex.path) bam_dir/$(inputs.sampleBamIndex.basename)
                mkdir variants_dir
                ln -s $(inputs.variantsForContamination.path) variants_dir/$(inputs.variantsForContamination.basename)
                ln -s $(inputs.variantsForContamination.path) variants_dir/$(inputs.variantsForContaminationIndex.basename)
                mkdir sites_dir
                ln -s $(inputs.sitesForContamination.path) sites_dir/$(inputs.sitesForContamination.basename)
                ln -s $(inputs.sitesForContaminationIndex.path) sites_dir/$(inputs.sitesForContaminationIndex.basename)
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                GetPileupSummaries \
                -I bam_dir/$(inputs.sampleBam.basename) \
                -V variants_dir/$(inputs.variantsForContamination.basename) \
                -L sites_dir/$(inputs.sitesForContamination.basename) \
                -O -pileups.table $(inputs.outputPrefix)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(120 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: HaplotypeCaller
    inputs:
      - id: inputBams
        type:
            items: File
            type: array
      - id: inputBamsIndex
        type:
            items: File
            type: array
      - id: outputPath
        type: string
      - id: referenceFasta
        type: File
      - id: referenceFastaIndex
        type: File
      - id: referenceFastaDict
        type: File
      - id: gvcf
        default: false
        type: boolean
      - id: emitRefConfidence
        type: string
      - id: dontUseSoftClippedBases
        default: false
        type: boolean
      - id: intervalList
        type:
          - items: File
            type: array
          - 'null'
      - id: excludeIntervalList
        type:
          - items: File
            type: array
          - 'null'
      - id: contamination
        type:
          - float
          - 'null'
      - id: dbsnpVCF
        type:
          - File
          - 'null'
      - id: dbsnpVCFIndex
        type:
          - File
          - 'null'
      - id: pedigree
        type:
          - File
          - 'null'
      - id: ploidy
        type:
          - int
          - 'null'
      - id: outputMode
        type:
          - string
          - 'null'
      - id: standardMinConfidenceThresholdForCalling
        type:
          - float
          - 'null'
      - id: javaXmxMb
        default: 4096
        type: int
      - id: memoryMb
        type:
          - int
          - 'null'
      - id: timeMinutes
        default: 400
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: outputVCF
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVCFIndex
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPath))"
                mkdir wd
                for FILE in $(inputs.inputBams.map(function(el) {return el.path}).join(" ")); do ln -s $FILE wd/\$(inputBams $FILE) ; done
                for FILE in $(inputs.inputBamsIndex.map(function(el) {return el.path}).join(" ")); do ln -s $FILE wd/\$(inputBamsIndex $FILE) ; done
                mkdir reference_dir
                ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
                ln -s $(inputs.referenceFastaDict.path) reference_dir/$(inputs.referenceFastaDict.basename)
                ln -s $(inputs.referenceFastaIndex.path) reference_dir/$(inputs.referenceFastaIndex.basename)
                gatk --java-options '-Xmx$(inputs.javaXmxMb)M -XX:ParallelGCThreads=1' \
                HaplotypeCaller \
                -R reference_dir/$(inputs.referenceFasta.basename) \
                -O $(inputs.outputPath) \
                (for FILE in $(inputs.inputBams.map(function(el) {return el.path}).join(" ")); do echo -- "-I wd/"\$(basename $FILE); done) \
                $(inputs.ploidy === null ? "" : "--sample-ploidy " + inputs.ploidy) \
                $(inputs.intervalList === null ? "" : "-L") $(inputs.intervalList.map(function(el) {return el.path}).join(" -L ")) \
                $(inputs.excludeIntervalList === null ? "" : "-XL") $(inputs.excludeIntervalList.map(function(el) {return el.path}).join(" -XL ")) \
                $(inputs.dbsnpVCF === null ? "" : "-D " + inputs.dbsnpVCF.path) \
                $(inputs.pedigree === null ? "" : "--pedigree " + inputs.pedigree.path) \
                $(inputs.contamination === null ? "" : "--contamination-fraction-per-sample-file " + inputs.contamination) \
                $(inputs.outputMode === null ? "" : "--output-mode " + inputs.outputMode) \
                --emit-ref-confidence $(inputs.emitRefConfidence) \
                $(inputs.dontUseSoftClippedBases ? "--dont-use-soft-clipped-bases" : "") \
                $(inputs.standardMinConfidenceThresholdForCalling === null ? "" : "--standard-min-confidence-threshold-for-calling " + inputs.standardMinConfidenceThresholdForCalling)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = "M";
            var value = parseInt(`${inputs.javaXmxMb + 512}`.match(/[0-9]+/g));
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(400 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: LearnReadOrientationModel
    inputs:
      - id: f1r2TarGz
        type:
            items: File
            type: array
      - id: javaXmx
        default: 12G
        type: string
      - id: memory
        default: 13G
        type: string
      - id: timeMinutes
        default: 120
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: artifactPriorsTable
        type: File
        outputBinding:
            glob: artifact-priors.tar.gz
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                LearnReadOrientationModel \
                -I $(inputs.f1r2TarGz.map(function(el) {return el.path}).join(" -I ")) \
                -O "artifact-priors.tar.gz"
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(120 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: MergeStats
    inputs:
      - id: stats
        type:
            items: File
            type: array
      - id: javaXmx
        default: 14G
        type: string
      - id: memory
        default: 15G
        type: string
      - id: timeMinutes
        default: 30
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: mergedStats
        type: File
        outputBinding:
            glob: merged.stats
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                MergeMutectStats \
                -stats $(inputs.stats.map(function(el) {return el.path}).join(" -stats ")) \
                -O "merged.stats"
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(30 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: ModelSegments
    inputs:
      - id: outputDir
        default: .
        type: string
      - id: outputPrefix
        type: string
      - id: denoisedCopyRatios
        type: File
      - id: allelicCounts
        type: File
      - id: minimumTotalAlleleCountCase
        type: int
      - id: maximumNumberOfSmoothingIterations
        default: 10
        type: int
      - id: normalAllelicCounts
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 10G
        type: string
      - id: memory
        default: 11G
        type: string
      - id: timeMinutes
        default: 60
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: hetrozygousAllelicCounts
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".hets.tsv")
      - id: copyRatioSegments
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".cr.seg")
      - id: copyRatioCBS
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".cr.igv.seg")
      - id: alleleFractionCBS
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".af.igv.seg")
      - id: unsmoothedModeledSegments
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelBegin.seg")
      - id: unsmoothedCopyRatioParameters
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelBegin.cr.param")
      - id: unsmoothedAlleleFractionParameters
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelBegin.af.param")
      - id: modeledSegments
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelFinal.seg")
      - id: copyRatioParameters
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelFinal.cr.param")
      - id: alleleFractionParameters
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelFinal.af.param")
      - id: normalHetrozygousAllelicCounts
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".hets.normal.tsv")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p $(inputs.outputDir)
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                ModelSegments \
                --denoised-copy-ratios $(inputs.denoisedCopyRatios.path) \
                --allelic-counts $(inputs.allelicCounts.path) \
                $(inputs.normalAllelicCounts === null ? "" : "--normal-allelic-counts " + inputs.normalAllelicCounts.path) \
                --minimum-total-allele-count-case $(inputs.minimumTotalAlleleCountCase) \
                --maximum-number-of-smoothing-iterations $(inputs.maximumNumberOfSmoothingIterations) \
                --output $(inputs.outputDir) \
                --output-prefix $(inputs.outputPrefix)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: MuTect2
    inputs:
      - id: inputBams
        type:
            items: File
            type: array
      - id: inputBamsIndex
        type:
            items: File
            type: array
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: outputVcf
        type: string
      - id: tumorSample
        type: string
      - id: f1r2TarGz
        default: f1r2.tar.gz
        type: string
      - id: intervals
        type:
            items: File
            type: array
      - id: outputStats
        type:
          - string
          - 'null'
      - id: normalSample
        type:
          - string
          - 'null'
      - id: germlineResource
        type:
          - File
          - 'null'
      - id: germlineResourceIndex
        type:
          - File
          - 'null'
      - id: panelOfNormals
        type:
          - File
          - 'null'
      - id: panelOfNormalsIndex
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 4G
        type: string
      - id: memory
        default: 5G
        type: string
      - id: timeMinutes
        default: 240
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: vcfFile
        type: File
        outputBinding:
            glob: $(inputs.outputVcf)
      - id: vcfFileIndex
        type: File
        outputBinding:
            glob: $(inputs.outputVcf + ".tbi")
      - id: f1r2File
        type: File
        outputBinding:
            glob: $(inputs.f1r2TarGz)
      - id: stats
        type: File
        outputBinding:
            glob: $(inputs.outputVcf + ".stats")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputVcf))"
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                Mutect2 \
                -R $(inputs.referenceFasta.path) \
                -I $(inputs.inputBams.map(function(el) {return el.path}).join(" -I ")) \
                -tumor $(inputs.tumorSample) \
                $(inputs.normalSample === null ? "" : "-normal " + inputs.normalSample) \
                $(inputs.germlineResource === null ? "" : "--germline-resource " + inputs.germlineResource.path) \
                $(inputs.panelOfNormals === null ? "" : "--panel-of-normals " + inputs.panelOfNormals.path) \
                --f1r2-tar-gz  $(inputs.f1r2TarGz) \
                -O $(inputs.outputVcf) \
                -L $(inputs.intervals.map(function(el) {return el.path}).join(" -L "))
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(240 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: PlotDenoisedCopyRatios
    inputs:
      - id: referenceFastaDict
        type: File
      - id: outputDir
        default: .
        type: string
      - id: outputPrefix
        type: string
      - id: standardizedCopyRatios
        type: File
      - id: denoisedCopyRatios
        type: File
      - id: minimumContigLength
        type:
          - int
          - 'null'
      - id: javaXmx
        default: 3G
        type: string
      - id: memory
        default: 4G
        type: string
      - id: timeMinutes
        default: 2
        type: int
      - id: dockerImage
        default: broadinstitute/gatk:4.1.8.0
        type: string
    outputs:
      - id: denoisedCopyRatiosPlot
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".denoised.png")
      - id: standardizedMedianAbsoluteDeviation
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".standardizedMAD.txt")
      - id: denoisedMedianAbsoluteDeviation
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".denoisedMAD.txt")
      - id: deltaMedianAbsoluteDeviation
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".deltaMAD.txt")
      - id: deltaScaledMedianAbsoluteDeviation
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".scaledDeltaMAD.txt")
      - id: denoisedCopyRatiosLimitedPlot
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".denoisedLimit4.png")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p $(inputs.outputDir)
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                PlotDenoisedCopyRatios \
                --standardized-copy-ratios $(inputs.standardizedCopyRatios.path) \
                --denoised-copy-ratios $(inputs.denoisedCopyRatios.path) \
                --sequence-dictionary $(inputs.referenceFastaDict.path) \
                $(inputs.minimumContigLength === null ? "" : "--minimum-contig-length " + inputs.minimumContigLength) \
                --output $(inputs.outputDir) \
                --output-prefix $(inputs.outputPrefix)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: broadinstitute/gatk:4.1.8.0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(2 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: PlotModeledSegments
    inputs:
      - id: referenceFastaDict
        type: File
      - id: outputDir
        default: .
        type: string
      - id: outputPrefix
        type: string
      - id: denoisedCopyRatios
        type: File
      - id: segments
        type: File
      - id: allelicCounts
        type: File
      - id: minimumContigLength
        type:
          - int
          - 'null'
      - id: javaXmx
        default: 3G
        type: string
      - id: memory
        default: 4G
        type: string
      - id: timeMinutes
        default: 2
        type: int
      - id: dockerImage
        default: broadinstitute/gatk:4.1.8.0
        type: string
    outputs:
      - id: modeledSegmentsPlot
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modeled.png")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p $(inputs.outputDir)
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                PlotModeledSegments \
                --denoised-copy-ratios $(inputs.denoisedCopyRatios.path) \
                --allelic-counts $(inputs.allelicCounts.path) \
                --segments $(inputs.segments.path) \
                --sequence-dictionary $(inputs.referenceFastaDict.path) \
                $(inputs.minimumContigLength === null ? "" : "--minimum-contig-length " + inputs.minimumContigLength) \
                --output $(inputs.outputDir) \
                --output-prefix $(inputs.outputPrefix)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: broadinstitute/gatk:4.1.8.0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(2 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: PreprocessIntervals
    inputs:
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: outputIntervalList
        default: bins.interval_list
        type: string
      - id: binLength
        type: int
      - id: padding
        type: int
      - id: intervalMergingRule
        default: OVERLAPPING_ONLY
        type: string
      - id: intervals
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 3G
        type: string
      - id: memory
        default: 4G
        type: string
      - id: timeMinutes
        type:
          - int
          - 'null'
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: intervalList
        type: File
        outputBinding:
            glob: $(inputs.outputIntervalList)
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputIntervalList))"
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                PreprocessIntervals \
                -R $(inputs.referenceFasta.path) \
                --sequence-dictionary $(inputs.referenceFastaDict.path) \
                --bin-length $(inputs.binLength) \
                --padding $(inputs.padding) \
                $(inputs.intervals === null ? "" : "-L " + inputs.intervals.path) \
                --interval-merging-rule $(inputs.intervalMergingRule) \
                -O $(inputs.outputIntervalList)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.referenceFasta.path.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3 * 6)  * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: SelectVariants
    inputs:
      - id: inputVcf
        type: File
      - id: inputVcfIndex
        type: File
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: outputPath
        default: output.vcf.gz
        type: string
      - id: intervals
        default: []
        type:
            items: File
            type: array
      - id: selectTypeToInclude
        type:
          - string
          - 'null'
      - id: javaXmx
        default: 4G
        type: string
      - id: memory
        default: 5G
        type: string
      - id: timeMinutes
        default: 60
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: outputVcf
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVcfIndex
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPath))"
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                SelectVariants \
                -R $(inputs.referenceFasta.path) \
                -V $(inputs.inputVcf.path) \
                $(inputs.selectTypeToInclude === null ? "" : "--select-type-to-include " + inputs.selectTypeToInclude) \
                $(inputs.intervals.length > 0 ? "-L" : "") $(inputs.intervals.map(function(el) {return el.path}).join(" -L ")) \
                -O $(inputs.outputPath)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: SplitNCigarReads
    inputs:
      - id: inputBam
        type: File
      - id: inputBamIndex
        type: File
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: outputBam
        type: string
      - id: intervals
        default: []
        type:
            items: File
            type: array
      - id: javaXmx
        default: 4G
        type: string
      - id: memory
        default: 5G
        type: string
      - id: timeMinutes
        default: 120
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: bam
        type: File
        outputBinding:
            glob: $(inputs.outputBam)
      - id: bamIndex
        type: File
        outputBinding:
            glob: $(inputs.outputBam.replace("\.bam$", ".bai") )
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputBam))"
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                SplitNCigarReads \
                -I $(inputs.inputBam.path) \
                -R $(inputs.referenceFasta.path) \
                -O $(inputs.outputBam) \
                $(inputs.intervals.length > 0 ? "-L" : "") $(inputs.intervals.map(function(el) {return el.path}).join(" -L "))
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(120 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: VariantEval
    inputs:
      - id: evalVcfs
        type:
            items: File
            type: array
      - id: evalVcfsIndex
        type:
            items: File
            type: array
      - id: comparisonVcfs
        default: []
        type:
            items: File
            type: array
      - id: comparisonVcfsIndex
        default: []
        type:
            items: File
            type: array
      - id: intervals
        default: []
        type:
            items: File
            type: array
      - id: outputPath
        default: eval.table
        type: string
      - id: doNotUseAllStandardModules
        default: false
        type: boolean
      - id: doNotUseAllStandardStratifications
        default: false
        type: boolean
      - id: evalModules
        default: []
        type:
            items: string
            type: array
      - id: stratificationModules
        default: []
        type:
            items: string
            type: array
      - id: samples
        default: []
        type:
            items: string
            type: array
      - id: mergeEvals
        default: false
        type: boolean
      - id: referenceFasta
        type:
          - File
          - 'null'
      - id: referenceFastaDict
        type:
          - File
          - 'null'
      - id: referenceFastaFai
        type:
          - File
          - 'null'
      - id: dbsnpVCF
        type:
          - File
          - 'null'
      - id: dbsnpVCFIndex
        type:
          - File
          - 'null'
      - id: javaXmx
        default: 4G
        type: string
      - id: memory
        default: 5G
        type: string
      - id: timeMinutes
        type:
          - int
          - 'null'
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: table
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPath))"
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                VariantEval \
                --output $(inputs.outputPath) \
                $(inputs.evalVcfs.length > 0 ? "--eval" : "") $(inputs.evalVcfs.map(function(el) {return el.path}).join(" --eval ")) \
                $(inputs.comparisonVcfs.length > 0 ? "--comparison" : "") $(inputs.comparisonVcfs.map(function(el) {return el.path}).join(" --comparison ")) \
                $(inputs.referenceFasta === null ? "" : "-R " + inputs.referenceFasta.path) \
                $(inputs.dbsnpVCF === null ? "" : "--dbsnp " + inputs.dbsnpVCF.path) \
                $(inputs.intervals.length > 0 ? "-L" : "") $(inputs.intervals.map(function(el) {return el.path}).join(" -L ")) \
                $(inputs.samples.length > 0 ? "--sample" : "") $(inputs.samples.join(" --sample ")) \
                $(inputs.doNotUseAllStandardModules ? "--do-not-use-all-standard-modules" : "") \
                $(inputs.doNotUseAllStandardStratifications ? "--do-not-use-all-standard-stratifications" : "") \
                $(inputs.evalModules.length > 0 ? "-EV" : "") $(inputs.evalModules.join(" -EV ")) \
                $(inputs.stratificationModules.length > 0 ? "-ST" : "") $(inputs.stratificationModules.join(" -ST ")) \
                $(inputs.mergeEvals ? "--merge-evals" : "")
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
      - class: ResourceRequirement
        coresMin: 1
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: '$(Math.ceil((function(size_of=0){(function () {var new_array =
            []; [ inputs.evalVcfs, inputs.comparisonVcfs, [inputs.referenceFasta ===
            null ? "" : inputs.referenceFasta.path, inputs.dbsnpVCF === null ? ""
            : inputs.dbsnpVCF.path].filter(function(element) { return element !==
            null })  ].forEach(function(value, index, obj) {value.forEach(function(sub_value,
            sub_index, sub_obj) {new_array.push(sub_value);});}); return new_array;})().forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3 * 20)  * 60)'
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: VariantFiltration
    inputs:
      - id: inputVcf
        type: File
      - id: inputVcfIndex
        type: File
      - id: referenceFasta
        type: File
      - id: referenceFastaDict
        type: File
      - id: referenceFastaFai
        type: File
      - id: outputPath
        default: filtered.vcf.gz
        type: string
      - id: filterArguments
        type:
            items: string
            type: array
      - id: intervals
        default: []
        type:
            items: File
            type: array
      - id: javaXmx
        default: 4G
        type: string
      - id: memory
        default: 5G
        type: string
      - id: timeMinutes
        default: 120
        type: int
      - id: dockerImage
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    outputs:
      - id: filteredVcf
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: filteredVcfIndex
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPath))"
                gatk --java-options '-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1' \
                VariantFiltration \
                -I $(inputs.inputVcf.path) \
                -R $(inputs.referenceFasta.path) \
                -O $(inputs.outputPath) \
                $(inputs.filterArguments.join(" ")) \
                $(inputs.intervals.length > 0 ? "-L" : "") $(inputs.intervals.map(function(el) {return el.path}).join(" -L "))
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(120 * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
