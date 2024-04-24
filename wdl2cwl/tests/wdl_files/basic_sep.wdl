version 1.1

workflow sepWorkflow {

    input {
        String delimiter
        Array[String] to_sep
    }
    output {
        String str_output = sep(delimiter, to_sep)
    }
}
