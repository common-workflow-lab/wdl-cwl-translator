version 1.0
task echo {
    input {
        Array[String] a_s
        Array[String] a_s2
    }

    command {
    }

    output {
        Array[String]+ out_s = flatten([a_s, a_s2])
    }
}

