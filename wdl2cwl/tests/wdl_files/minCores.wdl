version 1.0

task Echo {
    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
        cpu: "8"
    }
    command {
        echo "Hello world"
    }
    output {
        String out = stdout()
    }
}
