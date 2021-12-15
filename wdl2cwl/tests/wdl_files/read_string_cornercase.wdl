version 1.0


task ReadStringCornercase {
        input {
        File bam
        String expected_checksum
    }

    Int required_disk = ceil(size(bam, "G") * 1.1)

    command <<<
        cacheInvalidationRandomString=4

        echo Starting checksum generation...

        # calculate hash for alignment positions only (a reduced bam hash)
        calculated_checksum=$( samtools view -F 256 "~{bam}" | cut -f 1-11 | md5sum | awk '{print $1}' )
        echo Reduced checksum generation complete

        if [ "$calculated_checksum" == "~{expected_checksum}" ]
        then
             echo Computed and expected bam hashes match \( "$calculated_checksum" \)
             printf PASS > result.txt
        else
             echo Computed \( "$calculated_checksum" \) and expected \( "~{expected_checksum}" \) bam file hashes do not match
             printf FAIL > result.txt
        fi
    >>>

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
        cpu: 1
        memory: "3.75 GB"
        #disks: "local-disk ${required_disk} HDD"
    }
    output {
        File result = "read_string.txt"
    }
}
