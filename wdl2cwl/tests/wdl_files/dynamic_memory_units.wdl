version 1.0

task complex_memory_amount_and_units {
    input {
        String memory_units
        String memory_amount
    }
    runtime {
        memory: "~{memory_amount} ~{memory_units}"
    }
    command {
        echo "Hello world"
    }
    output {
        File out = stdout()
    }
}

task complex_memory_units {
    input {
        String memory_units
    }
    runtime {
        memory: "10 ~{memory_units}"
    }
    command {
        echo "Hello world"
    }
    output {
        File out = stdout()
    }
}

task complex_memory_amount {
    input {
        String memory_amount
    }
    runtime {
        memory: "~{memory_amount} KiB"
    }
    command {
        echo "Hello world"
    }
    output {
        File out = stdout()
    }
}


