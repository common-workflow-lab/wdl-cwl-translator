cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: AnnotateIntervals
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: referenceFasta
        doc: The reference fasta file.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: annotatedIntervalsPath
        doc: The location the output should be written to.
        default: intervals.annotated.tsv
        type: string
      - id: intervals
        doc: An interval list describinig the intervals to annotate.
        type: File
      - id: intervalMergingRule
        doc: Equivalent to gatk AnnotateIntervals' `--interval-merging-rule` option.
        default: OVERLAPPING_ONLY
        type: string
      - id: featureQueryLookahead
        doc: Equivalent to gatk AnnotateIntervals' `--feature-query-lookahead` option.
        default: 1000000
        type: int
      - id: mappabilityTrack
        doc: Equivalent to gatk AnnotateIntervals' `--mappability-track` option.
        type:
          - File
          - 'null'
      - id: segmentalDuplicationTrack
        doc: Equivalent to gatk AnnotateIntervals' `--segmenta-duplicarion-track`
            option.
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 2G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 3G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 5
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: annotatedIntervals
        doc: This is a tab-separated values (TSV) file with a SAM-style header containing
            a sequence dictionary, a row specifying the column headers for the contained
            annotations, and the corresponding entry rows.
        type: File
        outputBinding:
            glob: $(inputs.annotatedIntervalsPath)
  - cwlVersion: v1.2
    id: ApplyBQSR
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: inputBam
        doc: The BAM file which should be recalibrated.
        type: File
      - id: inputBamIndex
        doc: The input BAM file's index.
        type: File
      - id: outputBamPath
        doc: The location the resulting BAM file should be written.
        type: string
      - id: recalibrationReport
        doc: The BQSR report the be used for recalibration.
        type: File
      - id: sequenceGroupInterval
        doc: Bed files describing the regions to operate on.
        default: []
        type:
            items: File
            type: array
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: javaXmxMb
        doc: The maximum memory available to the program in megabytes. Should be lower
            than `memoryMb` to accommodate JVM overhead.
        default: 2048
        type: int
      - id: memoryMb
        doc: The amount of memory this job will use in megabytes.
        type:
          - int
          - 'null'
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: recalibratedBam
        doc: A BAM file containing the recalibrated read data.
        type: File
        outputBinding:
            glob: $(inputs.outputBamPath)
      - id: recalibratedBamIndex
        doc: Index of recalibrated BAM file.
        type: File
        outputBinding:
            glob: $(inputs.outputBamPath.replace("\\.bam$", ".bai") )
      - id: recalibratedBamMd5
        doc: MD5 of recalibrated BAM file.
        type: File
        outputBinding:
            glob: $(inputs.outputBamPath + ".md5")
  - cwlVersion: v1.2
    id: BaseRecalibrator
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: inputBam
        doc: The BAM file to generate a BQSR report for.
        type: File
      - id: inputBamIndex
        doc: The index of the input BAM file.
        type: File
      - id: recalibrationReportPath
        doc: The location to write the BQSR report to.
        type: string
      - id: sequenceGroupInterval
        doc: Bed files describing the regions to operate on.
        default: []
        type:
            items: File
            type: array
      - id: knownIndelsSitesVCFs
        doc: VCF files with known indels.
        default: []
        type:
            items: File
            type: array
      - id: knownIndelsSitesVCFIndexes
        doc: The indexed for the known variant VCFs.
        default: []
        type:
            items: File
            type: array
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: dbsnpVCF
        doc: A dbSNP VCF.
        type:
          - File
          - 'null'
      - id: dbsnpVCFIndex
        doc: The index for the dbSNP VCF.
        type:
          - File
          - 'null'
      - id: javaXmxMb
        doc: The maximum memory available to the program in megabytes. Should be lower
            than `memoryMb` to accommodate JVM overhead.
        default: 1024
        type: int
      - id: memoryMb
        doc: The amount of memory this job will use in megabytes.
        type:
          - int
          - 'null'
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: recalibrationReport
        doc: A GATK Report file with many tables.
        type: File
        outputBinding:
            glob: $(inputs.recalibrationReportPath)
  - cwlVersion: v1.2
    id: CalculateContamination
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: tumorPileups
        doc: The pileup summary of a tumor/case sample.
        type: File
      - id: normalPileups
        doc: The pileup summary of the normal/control sample.
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 12G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 13G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 180
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: contaminationTable
        doc: Table with fractions of reads from cross-sample contamination.
        type: File
        outputBinding:
            glob: contamination.table
      - id: mafTumorSegments
        doc: Tumor segments table.
        type: File
        outputBinding:
            glob: segments.table
  - cwlVersion: v1.2
    id: CallCopyRatioSegments
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: outputPrefix
        doc: The prefix for the output files.
        type: string
      - id: copyRatioSegments
        doc: The copy ratios file generated by gatk ModelSegments.
        type: File
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 2G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 3G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 2
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: calledSegments
        doc: This is a tab-separated values (TSV) file with a SAM-style header containing
            a read group sample name, a sequence dictionary, a row specifying the
            column headers contained in CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn,
            and the corresponding entry rows.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".called.seg")
      - id: calledSegmentsIgv
        doc: This is a tab-separated values (TSV) file with CBS-format column headers
            and the corresponding entry rows that can be plotted using IGV.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".called.igv.seg")
  - cwlVersion: v1.2
    id: CollectAllelicCounts
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: allelicCountsPath
        doc: The path the output should be written to.
        default: allelic_counts.tsv
        type: string
      - id: commonVariantSites
        doc: Interval list or vcf of common variant sites (to retrieve the allelic
            counts for).
        type: File
      - id: inputBam
        doc: The BAM file to generate counts for.
        type: File
      - id: inputBamIndex
        doc: The index of the input BAM file.
        type: File
      - id: referenceFasta
        doc: The reference fasta file.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: commonVariantSitesIndex
        doc: The index for commonVariantSites.
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 10G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 11G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: allelicCounts
        doc: This is a tab-separated values (TSV) file with a SAM-style header containing
            a read group sample name, a sequence dictionary, a row specifying the
            column headers contained in AllelicCountCollection.AllelicCountTableColumn,
            and the corresponding entry rows.
        type: File
        outputBinding:
            glob: $(inputs.allelicCountsPath)
  - cwlVersion: v1.2
    id: CollectReadCounts
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.inputBam.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3 * 5)  * 60)
    inputs:
      - id: countsPath
        doc: The location the output should be written to.
        default: readcounts.hdf5
        type: string
      - id: intervals
        doc: The intervals to collect counts for.
        type: File
      - id: inputBam
        doc: The BAM file to determine the coverage for.
        type: File
      - id: inputBamIndex
        doc: The input BAM file's index.
        type: File
      - id: referenceFasta
        doc: The reference fasta file.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: intervalMergingRule
        doc: Equivalent to gatk CollectReadCounts' `--interval-merging-rule` option.
        default: OVERLAPPING_ONLY
        type: string
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 7G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 8G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: counts
        doc: Read counts at specified intervals.
        type: File
        outputBinding:
            glob: $(inputs.countsPath)
  - cwlVersion: v1.2
    id: CombineGVCFs
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.gvcfFiles.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3 * 8)  * 60)
    inputs:
      - id: gvcfFiles
        doc: The GVCF files to be combined.
        type:
            items: File
            type: array
      - id: gvcfFilesIndex
        doc: The indexes for the GVCF files.
        type:
            items: File
            type: array
      - id: intervals
        doc: Bed files or interval lists describing the regions to operate on.
        default: []
        type:
            items: File
            type: array
      - id: outputPath
        doc: The location the combined GVCF should be written to.
        type: string
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 4G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 5G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.gvcfFiles.length == 0) {throw "gvcfFiles must contain
            at least one item.";} else { return "";}}
      - valueFrom: ${if (inputs.gvcfFilesIndex.length == 0) {throw "gvcfFilesIndex
            must contain at least one item.";} else { return "";}}
    outputs:
      - id: outputVcf
        doc: A combined multi-sample gVCF.
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVcfIndex
        doc: Index of the output file.
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
  - cwlVersion: v1.2
    id: CombineVariants
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(180 * 60)
    inputs:
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: genotypeMergeOption
        doc: Equivalent to CombineVariants' `--genotypemergeoption` option.
        default: UNIQUIFY
        type: string
      - id: filteredRecordsMergeType
        doc: Equivalent to CombineVariants' `--filteredrecordsmergetype` option.
        default: KEEP_IF_ANY_UNFILTERED
        type: string
      - id: identifiers
        doc: The sample identifiers in the same order as variantVcfs.
        type:
            items: string
            type: array
      - id: variantVcfs
        doc: The input VCF files in the same order as identifiers.
        type:
            items: File
            type: array
      - id: variantIndexes
        doc: The indexes of the input VCF files.
        type:
            items: File
            type: array
      - id: outputPath
        doc: The location the output should be written to.
        type: string
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 12G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 13G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 180
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: broadinstitute/gatk3:3.8-1
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.identifiers.length == 0) {throw "identifiers must
            contain at least one item.";} else { return "";}}
      - valueFrom: ${if (inputs.variantVcfs.length == 0) {throw "variantVcfs must
            contain at least one item.";} else { return "";}}
      - valueFrom: ${if (inputs.variantIndexes.length == 0) {throw "variantIndexes
            must contain at least one item.";} else { return "";}}
    outputs:
      - id: combinedVcf
        doc: Combined VCF file.
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: combinedVcfIndex
        doc: Index of combined VCF file.
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
  - cwlVersion: v1.2
    id: CreateReadCountPanelOfNormals
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(5 * 60)
    inputs:
      - id: PONpath
        doc: The location the PON should be written to.
        default: PON.hdf5
        type: string
      - id: readCountsFiles
        doc: The read counts files as generated by CollectReadCounts.
        type:
            items: File
            type: array
      - id: annotatedIntervals
        doc: An annotation set of intervals as generated by AnnotateIntervals. If
            provided, explicit GC correction will be performed.
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 7G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 8G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 5
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: broadinstitute/gatk:4.1.8.0
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.readCountsFiles.length == 0) {throw "readCountsFiles
            must contain at least one item.";} else { return "";}}
    outputs:
      - id: PON
        doc: Panel-of-normals file.
        type: File
        outputBinding:
            glob: $(inputs.PONpath)
  - cwlVersion: v1.2
    id: DenoiseReadCounts
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(5 * 60)
    inputs:
      - id: readCounts
        doc: The read counts file as generated by CollectReadCounts.
        type: File
      - id: outputPrefix
        doc: The prefix for the output files.
        type: string
      - id: PON
        doc: A panel of normals as generated by CreateReadCountPanelOfNormals.
        type:
          - File
          - 'null'
      - id: annotatedIntervals
        doc: An annotated set of intervals as generated by AnnotateIntervals. Will
            be ignored if PON is provided.
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 4G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 5G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 5
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: standardizedCopyRatios
        doc: This is a tab-separated values (TSV) file with a SAM-style header containing
            a read group sample name, a sequence dictionary, a row specifying the
            column headers contained in CopyRatioCollection.CopyRatioTableColumn,
            and the corresponding entry rows.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".standardizedCR.tsv")
      - id: denoisedCopyRatios
        doc: This is a tab-separated values (TSV) file with a SAM-style header containing
            a read group sample name, a sequence dictionary, a row specifying the
            column headers contained in CopyRatioCollection.CopyRatioTableColumn,
            and the corresponding entry rows.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".denoisedCR.tsv")
  - cwlVersion: v1.2
    id: FilterMutectCalls
    class: CommandLineTool
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
                $("--unique-alt-read-count " + inputs.uniqueAltReadCount) \
                $("-stats " + inputs.mutect2Stats.path) \
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(60 * 60)
    inputs:
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: unfilteredVcf
        doc: An unfiltered VCF file as produced by Mutect2.
        type: File
      - id: unfilteredVcfIndex
        doc: The index of the unfiltered VCF file.
        type: File
      - id: outputVcf
        doc: The location the filtered VCF file should be written.
        type: string
      - id: uniqueAltReadCount
        doc: Equivalent to FilterMutectCalls' `--unique-alt-read-count` option.
        default: 4
        type: int
      - id: mutect2Stats
        doc: Equivalent to FilterMutectCalls' `-stats` option.
        type: File
      - id: contaminationTable
        doc: Equivalent to FilterMutectCalls' `--contamination-table` option.
        type:
          - File
          - 'null'
      - id: mafTumorSegments
        doc: Equivalent to FilterMutectCalls' `--tumor-segmentation` option.
        type:
          - File
          - 'null'
      - id: artifactPriors
        doc: Equivalent to FilterMutectCalls' `--ob-priors` option.
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 12G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 13G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 60
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: filteredVcf
        doc: VCF file with filtered variants from a Mutect2 VCF callset.
        type: File
        outputBinding:
            glob: $(inputs.outputVcf)
      - id: filteredVcfIndex
        doc: Index of output VCF file.
        type: File
        outputBinding:
            glob: $(inputs.outputVcf + ".tbi")
      - id: filteringStats
        doc: The output filtering stats file.
        type: File
        outputBinding:
            glob: filtering.stats
  - cwlVersion: v1.2
    id: GatherBqsrReports
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(1 * 60)
    inputs:
      - id: inputBQSRreports
        doc: The BQSR reports to be merged.
        type:
            items: File
            type: array
      - id: outputReportPath
        doc: The location of the combined BQSR report.
        type: string
      - id: javaXmxMb
        doc: The maximum memory available to the program in megabytes. Should be lower
            than `memory` to accommodate JVM overhead.
        default: 256
        type: int
      - id: memoryMb
        doc: The amount of memory this job will use in megabytes.
        type:
          - int
          - 'null'
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 1
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: outputBQSRreport
        doc: Single file with scattered BQSR recalibration reports gathered into one.
        type: File
        outputBinding:
            glob: $(inputs.outputReportPath)
  - cwlVersion: v1.2
    id: GenomicsDBImport
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(180 * 60)
    inputs:
      - id: gvcfFiles
        doc: The gvcfFiles to be merged.
        type:
            items: File
            type: array
      - id: gvcfFilesIndex
        doc: Indexes for the gvcfFiles.
        type:
            items: File
            type: array
      - id: intervals
        doc: intervals over which to operate.
        type:
            items: File
            type: array
      - id: genomicsDBWorkspacePath
        doc: Where the genomicsDB files should be stored.
        default: genomics_db
        type: string
      - id: genomicsDBTarFile
        doc: Where the .tar file containing the genomicsDB should be stored.
        default: genomics_db.tar.gz
        type: string
      - id: tmpDir
        doc: Alternate temporary directory in case there is not enough space. Must
            be mounted when using containers.
        type:
          - string
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 4G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 5G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 180
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.intervals.length == 0) {throw "intervals must contain
            at least one item.";} else { return "";}}
    outputs:
      - id: genomicsDbTarArchive
        doc: Imported VCFs to GenomicsDB file.
        type: File
        outputBinding:
            glob: $(inputs.genomicsDBTarFile)
  - cwlVersion: v1.2
    id: GenotypeGVCFs
    class: CommandLineTool
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
                $(inputs.intervals !== null ? "--only-output-calls-starting-in-intervals" : "") \
                $(inputs.intervals !== null ? "-L" : "") $(inputs.intervals.map(function(el) {return el.path}).join(" -L "))
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(120 * 60)
    inputs:
      - id: gvcfFile
        doc: The GVCF file to be genotyped.
        type: File
      - id: gvcfFileIndex
        doc: The index of the input GVCF file.
        type: File
      - id: outputPath
        doc: The location to write the output VCF file to.
        type: string
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: annotationGroups
        doc: Which annotation groups will be used for the annotation.
        default:
          - StandardAnnotation
        type:
            items: string
            type: array
      - id: intervals
        doc: Bed files or interval lists describing the regions to operate on.
        type:
          - items: File
            type: array
          - 'null'
      - id: dbsnpVCF
        doc: A dbSNP VCF.
        type:
          - File
          - 'null'
      - id: dbsnpVCFIndex
        doc: The index for the dbSNP VCF.
        type:
          - File
          - 'null'
      - id: pedigree
        doc: Pedigree file for determining the population "founders".
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 6G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 7G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: outputVCF
        doc: 'A final VCF in which all samples have been jointly genotyped. '
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVCFIndex
        doc: Index of final VCF file.
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
  - cwlVersion: v1.2
    id: GetPileupSummaries
    class: CommandLineTool
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
                -O $("-pileups.table" + inputs.outputPrefix)
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(120 * 60)
    inputs:
      - id: sampleBam
        doc: A BAM file for which a pileup should be created.
        type: File
      - id: sampleBamIndex
        doc: The index of the input BAM file.
        type: File
      - id: variantsForContamination
        doc: A VCF file with common variants.
        type: File
      - id: variantsForContaminationIndex
        doc: The index for the common variants VCF file.
        type: File
      - id: sitesForContamination
        doc: A bed file describing regions to operate on.
        type: File
      - id: sitesForContaminationIndex
        doc: The index for the bed file.
        type: File
      - id: outputPrefix
        doc: The prefix for the ouput.
        type: string
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 12G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 13G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: pileups
        doc: Pileup metrics for inferring contamination.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "-pileups.table")
  - cwlVersion: v1.2
    id: HaplotypeCaller
    class: CommandLineTool
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
                $(inputs.intervalList !== null ? "-L" : "") $(inputs.intervalList.map(function(el) {return el.path}).join(" -L ")) \
                $(inputs.excludeIntervalList !== null ? "-XL" : "") $(inputs.excludeIntervalList.map(function(el) {return el.path}).join(" -XL ")) \
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(400 * 60)
    inputs:
      - id: inputBams
        doc: The BAM files on which to perform variant calling.
        type:
            items: File
            type: array
      - id: inputBamsIndex
        doc: The indexes for the input BAM files.
        type:
            items: File
            type: array
      - id: outputPath
        doc: The location to write the output to.
        type: string
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaIndex
        doc: The index for the reference fasta file.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: gvcf
        doc: Whether the output should be a gvcf.
        default: false
        type: boolean
      - id: emitRefConfidence
        doc: "Whether to include reference calls. Three modes: 'NONE', 'BP_RESOLUTION'\
            \ and 'GVCF'."
        type: string
      - id: dontUseSoftClippedBases
        doc: Do not use soft-clipped bases. Should be 'true' for RNA variant calling.
        default: false
        type: boolean
      - id: intervalList
        doc: Bed files or interval lists describing the regions to operate on.
        type:
          - items: File
            type: array
          - 'null'
      - id: excludeIntervalList
        doc: Bed files or interval lists describing the regions to NOT operate on.
        type:
          - items: File
            type: array
          - 'null'
      - id: contamination
        doc: Equivalent to HaplotypeCaller's `-contamination` option.
        type:
          - float
          - 'null'
      - id: dbsnpVCF
        doc: A dbSNP VCF.
        type:
          - File
          - 'null'
      - id: dbsnpVCFIndex
        doc: The index for the dbSNP VCF.
        type:
          - File
          - 'null'
      - id: pedigree
        doc: Pedigree file for determining the population "founders".
        type:
          - File
          - 'null'
      - id: ploidy
        doc: The ploidy with which the variants should be called.
        type:
          - int
          - 'null'
      - id: outputMode
        doc: Specifies which type of calls we should output. Same as HaplotypeCaller's
            `--output-mode` option.
        type:
          - string
          - 'null'
      - id: standardMinConfidenceThresholdForCalling
        doc: Confidence threshold used for calling variants.
        type:
          - float
          - 'null'
      - id: javaXmxMb
        doc: The maximum memory available to the program in megabytes. Should be lower
            than `memoryMb` to accommodate JVM overhead.
        default: 4096
        type: int
      - id: memoryMb
        doc: The amount of memory this job will use in megabytes.
        type:
          - int
          - 'null'
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 400
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.inputBams.length == 0) {throw "inputBams must contain
            at least one item.";} else { return "";}}
      - valueFrom: ${if (inputs.inputBamsIndex.length == 0) {throw "inputBamsIndex
            must contain at least one item.";} else { return "";}}
      - valueFrom: ${if (inputs.intervalList.length == 0) {throw "intervalList must
            contain at least one item.";} else { return "";}}
      - valueFrom: ${if (inputs.excludeIntervalList.length == 0) {throw "excludeIntervalList
            must contain at least one item.";} else { return "";}}
    outputs:
      - id: outputVCF
        doc: Raw, unfiltered SNP and indel calls.
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVCFIndex
        doc: Index of output VCF.
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
  - cwlVersion: v1.2
    id: LearnReadOrientationModel
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(120 * 60)
    inputs:
      - id: f1r2TarGz
        doc: A f1r2TarGz file outputed by mutect2.
        type:
            items: File
            type: array
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 12G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 13G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.f1r2TarGz.length == 0) {throw "f1r2TarGz must contain
            at least one item.";} else { return "";}}
    outputs:
      - id: artifactPriorsTable
        doc: Maximum likelihood estimates of artifact prior probabilities in the orientation
            bias mixture model filter.
        type: File
        outputBinding:
            glob: artifact-priors.tar.gz
  - cwlVersion: v1.2
    id: MergeStats
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(30 * 60)
    inputs:
      - id: stats
        doc: Statistics files to be merged.
        type:
            items: File
            type: array
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 14G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 15G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 30
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.stats.length == 0) {throw "stats must contain at least
            one item.";} else { return "";}}
    outputs:
      - id: mergedStats
        doc: Merged stats from scattered Mutect2 runs.
        type: File
        outputBinding:
            glob: merged.stats
  - cwlVersion: v1.2
    id: ModelSegments
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
    inputs:
      - id: outputDir
        doc: The directory to write the ouput to.
        default: .
        type: string
      - id: outputPrefix
        doc: The prefix of the output files. Should not include directories.
        type: string
      - id: denoisedCopyRatios
        doc: The denoised copy ratios as generated by DenoiseReadCounts.
        type: File
      - id: allelicCounts
        doc: The allelicCounts as generate by CollectAllelicCounts.
        type: File
      - id: minimumTotalAlleleCountCase
        doc: Equivalent to gatk ModelSeqments' `--minimum-total-allele-count-case`
            option.
        type: int
      - id: maximumNumberOfSmoothingIterations
        doc: Equivalent to gatk ModelSeqments' `--maximum-number-of-smoothing-iterations`
            option.
        default: 10
        type: int
      - id: normalAllelicCounts
        doc: The allelicCounts as generate by CollectAllelicCounts for a matched normal.
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 10G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 11G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 60
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: hetrozygousAllelicCounts
        doc: Allelic-counts file containing the counts at sites genotyped as heterozygous
            in the case sample.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".hets.tsv")
      - id: copyRatioSegments
        doc: It contains the segments from the .modelFinal.seg file converted to a
            format suitable for input to CallCopyRatioSegments.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".cr.seg")
      - id: copyRatioCBS
        doc: The posterior medians of the log2 copy ratio.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".cr.igv.seg")
      - id: alleleFractionCBS
        doc: Minor-allele fraction.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".af.igv.seg")
      - id: unsmoothedModeledSegments
        doc: The initial modeled-segments result before segmentation smoothing.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelBegin.seg")
      - id: unsmoothedCopyRatioParameters
        doc: The initial copy-ratio-model global-parameter result before segmentation
            smoothing.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelBegin.cr.param")
      - id: unsmoothedAlleleFractionParameters
        doc: The initial allele-fraction-model global-parameter result before segmentation
            smoothing.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelBegin.af.param")
      - id: modeledSegments
        doc: The final modeled-segments result after segmentation smoothing.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelFinal.seg")
      - id: copyRatioParameters
        doc: The final copy-ratio-model global-parameter result after segmentation
            smoothing.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelFinal.cr.param")
      - id: alleleFractionParameters
        doc: The final allele-fraction-model global-parameter result after segmentation
            smoothing.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modelFinal.af.param")
      - id: normalHetrozygousAllelicCounts
        doc: Allelic-counts file containing the counts at sites genotyped as heterozygous
            in the matched-normal sample.
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".hets.normal.tsv")
  - cwlVersion: v1.2
    id: MuTect2
    class: CommandLineTool
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
                $("--f1r2-tar-gz " + inputs.f1r2TarGz) \
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(240 * 60)
    inputs:
      - id: inputBams
        doc: The BAM files on which to perform variant calling.
        type:
            items: File
            type: array
      - id: inputBamsIndex
        doc: The indexes for the input BAM files.
        type:
            items: File
            type: array
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: outputVcf
        doc: The location to write the output VCF file to.
        type: string
      - id: tumorSample
        doc: The name of the tumor/case sample.
        type: string
      - id: f1r2TarGz
        doc: Equivalent to Mutect2's `--f1r2-tar-gz` option.
        default: f1r2.tar.gz
        type: string
      - id: intervals
        doc: Bed files describing the regiosn to operate on.
        type:
            items: File
            type: array
      - id: outputStats
        doc: The location the output statistics should be written to.
        type:
          - string
          - 'null'
      - id: normalSample
        doc: The name of the normal/control sample.
        type:
          - string
          - 'null'
      - id: germlineResource
        doc: Equivalent to Mutect2's `--germline-resource` option.
        type:
          - File
          - 'null'
      - id: germlineResourceIndex
        doc: The index for the germline resource.
        type:
          - File
          - 'null'
      - id: panelOfNormals
        doc: Equivalent to Mutect2's `--panel-of-normals` option.
        type:
          - File
          - 'null'
      - id: panelOfNormalsIndex
        doc: The index for the panel of normals.
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 4G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 5G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 240
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.inputBams.length == 0) {throw "inputBams must contain
            at least one item.";} else { return "";}}
      - valueFrom: ${if (inputs.inputBamsIndex.length == 0) {throw "inputBamsIndex
            must contain at least one item.";} else { return "";}}
      - valueFrom: ${if (inputs.intervals.length == 0) {throw "intervals must contain
            at least one item.";} else { return "";}}
    outputs:
      - id: vcfFile
        doc: Somatic SNVs and indels called via local assembly of haplotypes.
        type: File
        outputBinding:
            glob: $(inputs.outputVcf)
      - id: vcfFileIndex
        doc: Index for Mutect2 VCF.
        type: File
        outputBinding:
            glob: $(inputs.outputVcf + ".tbi")
      - id: f1r2File
        doc: Contains information that can then be passed to LearnReadOrientationModel,
            which generate an artifact prior table for each tumor sample for FilterMutectCalls
            to use.
        type: File
        outputBinding:
            glob: $(inputs.f1r2TarGz)
      - id: stats
        doc: Stats file.
        type: File
        outputBinding:
            glob: $(inputs.outputVcf + ".stats")
  - cwlVersion: v1.2
    id: PlotDenoisedCopyRatios
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(2 * 60)
    inputs:
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file used
            for the analyses.
        type: File
      - id: outputDir
        doc: The directory to write the ouput to.
        default: .
        type: string
      - id: outputPrefix
        doc: The prefix of the output files. Should not include directories.
        type: string
      - id: standardizedCopyRatios
        doc: The standardized copy ratios as generated by DenoiseReadCounts.
        type: File
      - id: denoisedCopyRatios
        doc: The denoised copy ratios as generated by DenoiseReadCounts.
        type: File
      - id: minimumContigLength
        doc: The minimum length for a contig to be included in the plots.
        type:
          - int
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 3G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 2
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: broadinstitute/gatk:4.1.8.0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: denoisedCopyRatiosPlot
        doc: Plot showing the entire range of standardized and denoised copy ratios.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".denoised.png")
      - id: standardizedMedianAbsoluteDeviation
        doc: Standardized median absolute deviation copy ratios.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".standardizedMAD.txt")
      - id: denoisedMedianAbsoluteDeviation
        doc: Denoised median absolute deviation copy ratios.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".denoisedMAD.txt")
      - id: deltaMedianAbsoluteDeviation
        doc: The change between `standardizedMedianAbsoluteDeviation` & `denoisedMedianAbsoluteDeviation`.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".deltaMAD.txt")
      - id: deltaScaledMedianAbsoluteDeviation
        doc: The change between `standardizedMedianAbsoluteDeviation` & `denoisedMedianAbsoluteDeviation`
            scaled by standardized MAD.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".scaledDeltaMAD.txt")
      - id: denoisedCopyRatiosLimitedPlot
        doc: Plot showing the standardized and denoised copy ratios limited to ratios
            within [0, 4].
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".denoisedLimit4.png")
  - cwlVersion: v1.2
    id: PlotModeledSegments
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(2 * 60)
    inputs:
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file used
            for the analyses.
        type: File
      - id: outputDir
        doc: The directory to write the ouput to.
        default: .
        type: string
      - id: outputPrefix
        doc: The prefix of the output files. Should not include directories.
        type: string
      - id: denoisedCopyRatios
        doc: The denoised copy ratios as generated by DenoiseReadCounts.
        type: File
      - id: segments
        doc: The modeled segments as generated by ModelSegments.
        type: File
      - id: allelicCounts
        doc: The hetrozygous allelic counts as generated by ModelSegments.
        type: File
      - id: minimumContigLength
        doc: The minimum length for a contig to be included in the plots.
        type:
          - int
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 3G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 2
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: broadinstitute/gatk:4.1.8.0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: modeledSegmentsPlot
        doc: This plot shows the input denoised copy ratios and/or alternate-allele
            fractions as points, as well as box plots for the available posteriors
            in each segment.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/" + inputs.outputPrefix + ".modeled.png")
  - cwlVersion: v1.2
    id: PreprocessIntervals
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.referenceFasta.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3 * 6)  * 60)
    inputs:
      - id: referenceFasta
        doc: The reference fasta file.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: outputIntervalList
        doc: The location the output should be written to.
        default: bins.interval_list
        type: string
      - id: binLength
        doc: The size of the bins to be created. Should be 0 for targeted/exome sequencing.
        type: int
      - id: padding
        doc: The padding to be added to the bins. Should be 0 if contiguos binning
            is used, eg with WGS.
        type: int
      - id: intervalMergingRule
        doc: Equivalent to gatk PreprocessIntervals' `--interval-merging-rule` option.
        default: OVERLAPPING_ONLY
        type: string
      - id: intervals
        doc: Bed files describing the regiosn to operate on.
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 3G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: intervalList
        doc: Preprocessed Picard interval-list file.
        type: File
        outputBinding:
            glob: $(inputs.outputIntervalList)
  - cwlVersion: v1.2
    id: SelectVariants
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
    inputs:
      - id: inputVcf
        doc: The VCF input file.
        type: File
      - id: inputVcfIndex
        doc: The input VCF file's index.
        type: File
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: outputPath
        doc: The location the output VCF file should be written.
        default: output.vcf.gz
        type: string
      - id: intervals
        doc: Bed files or interval lists describing the regions to operate on.
        default: []
        type:
            items: File
            type: array
      - id: selectTypeToInclude
        doc: Select only a certain type of variants from the input file.
        type:
          - string
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 4G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 5G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 60
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: outputVcf
        doc: A new VCF file containing the selected subset of variants.
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVcfIndex
        doc: Index of the new output VCF file.
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
  - cwlVersion: v1.2
    id: SplitNCigarReads
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(120 * 60)
    inputs:
      - id: inputBam
        doc: The BAM file for which spliced reads should be split.
        type: File
      - id: inputBamIndex
        doc: The input BAM file's index.
        type: File
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: outputBam
        doc: The location the output BAM file should be written.
        type: string
      - id: intervals
        doc: Bed files or interval lists describing the regions to operate on.
        default: []
        type:
            items: File
            type: array
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 4G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 5G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: bam
        doc: BAM file with reads split at N CIGAR elements and CIGAR strings updated.
        type: File
        outputBinding:
            glob: $(inputs.outputBam)
      - id: bamIndex
        doc: Index of output BAM file.
        type: File
        outputBinding:
            glob: $(inputs.outputBam.replace("\\.bam$", ".bai") )
  - cwlVersion: v1.2
    id: VariantEval
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
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
    inputs:
      - id: evalVcfs
        doc: Variant sets to evaluate.
        type:
            items: File
            type: array
      - id: evalVcfsIndex
        doc: Indexes for the variant sets.
        type:
            items: File
            type: array
      - id: comparisonVcfs
        doc: Compare set vcfs.
        default: []
        type:
            items: File
            type: array
      - id: comparisonVcfsIndex
        doc: Indexes for the compare sets.
        default: []
        type:
            items: File
            type: array
      - id: intervals
        doc: Bed files or interval lists describing the regions to operate on.
        default: []
        type:
            items: File
            type: array
      - id: outputPath
        doc: The location the output table should be written.
        default: eval.table
        type: string
      - id: doNotUseAllStandardModules
        doc: Do not use the standard modules by default (instead, only those that
            are specified with the evalModules option).
        default: false
        type: boolean
      - id: doNotUseAllStandardStratifications
        doc: Do not use the standard stratification modules by default (instead, only
            those that are specified with the stratificationModules option).
        default: false
        type: boolean
      - id: evalModules
        doc: One or more specific eval modules to apply to the eval track(s) (in addition
            to the standard modules, unless doNotUseAllStandardModules=true).
        default: []
        type:
            items: string
            type: array
      - id: stratificationModules
        doc: One or more specific stratification modules to apply to the eval track(s)
            (in addition to the standard stratifications, unless doNotUseAllStandardStratifications=true).
        default: []
        type:
            items: string
            type: array
      - id: samples
        doc: Derive eval and comp contexts using only these sample genotypes, when
            genotypes are available in the original context.
        default: []
        type:
            items: string
            type: array
      - id: mergeEvals
        doc: If provided, all evalVcf tracks will be merged into a single eval track.
        default: false
        type: boolean
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type:
          - File
          - 'null'
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type:
          - File
          - 'null'
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type:
          - File
          - 'null'
      - id: dbsnpVCF
        doc: A dbSNP VCF.
        type:
          - File
          - 'null'
      - id: dbsnpVCFIndex
        doc: The index for the dbSNP VCF.
        type:
          - File
          - 'null'
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 4G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 5G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: table
        doc: Evaluation tables detailing the results of the eval modules which were
            applied.
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
  - cwlVersion: v1.2
    id: VariantFiltration
    class: CommandLineTool
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(120 * 60)
    inputs:
      - id: inputVcf
        doc: The VCF to be filtered.
        type: File
      - id: inputVcfIndex
        doc: The input VCF file's index.
        type: File
      - id: referenceFasta
        doc: The reference fasta file which was also used for mapping.
        type: File
      - id: referenceFastaDict
        doc: The sequence dictionary associated with the reference fasta file.
        type: File
      - id: referenceFastaFai
        doc: The index for the reference fasta file.
        type: File
      - id: outputPath
        doc: The location the output VCF file should be written.
        default: filtered.vcf.gz
        type: string
      - id: filterArguments
        doc: "Arguments that should be used for the filter. For example: ['--filter-name',\
            \ 'my_filter', '--filter-expression', 'AB<0.2']."
        type:
            items: string
            type: array
      - id: intervals
        doc: Bed files or interval lists describing the regions to operate on.
        default: []
        type:
            items: File
            type: array
      - id: javaXmx
        doc: The maximum memory available to the program. Should be lower than `memory`
            to accommodate JVM overhead.
        default: 4G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 5G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.filterArguments.length == 0) {throw "filterArguments
            must contain at least one item.";} else { return "";}}
    outputs:
      - id: filteredVcf
        doc: A filtered VCF in which passing variants are annotated as PASS and failing
            variants are annotated with the name(s) of the filter(s) they failed.
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: filteredVcfIndex
        doc: Index of filtered VCF.
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
