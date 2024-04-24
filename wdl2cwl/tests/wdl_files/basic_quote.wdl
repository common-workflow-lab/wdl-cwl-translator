# From https://github.com/DataBiosphere/wdl-conformance-tests/blob/c87b62b4f460e009fd42edec13669c4db14cf90c/tests/basic_quote/basic_quote.wdl

version 1.1

workflow quoteWorkflow {
    input {
        Array[String] str_arr
        Array[Int] int_arr
        Array[Float] float_arr
        Array[Boolean] bool_arr
        Array[File] file_arr
    }
    call file_quote {input: file_arr = file_arr}
    output {
        Array[String] str_output = quote(str_arr)
        Array[String] int_output = quote(int_arr)
        Array[String] float_output = quote(float_arr)
        Array[String] bool_output = quote(bool_arr)
        File file_output = file_quote.out
    }
}

task file_quote {
    input {
        Array[File] file_arr
    }

    command <<<
        echo ~{sep(' ', quote(file_arr))} >> output.txt
    >>>

    output {
        File out = "output.txt"
    }
}
