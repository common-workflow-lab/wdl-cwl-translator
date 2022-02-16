version 1.0
# mixing "default" with other placeholder options is not allowed in 1.1

task placeholder_options_test {
  Boolean b = true
  Array[Int] numbers = [1, 2, 3]

  input {
    String? absent
    Array[Int]? missing_numbers
  }

  command <<<
    echo ~{true="true single quote: 'foo' " false="error" b}
    echo ~{true='true double quote: "foo" ' false="error" b}
    echo ~{true='true mixed quotes: "foo" \'bar\' ' false="error" b}
    echo ~{true="true mixed quotes: \"foo\" 'bar' " false="error" b}
    echo ~{default="default single 'quote'" absent}
    echo ~{default='default double "quote"' absent}
    echo ~{default='default "mixed" \'quotes\'' absent}
    echo ~{default="default \"mixed\" 'quotes'" absent}
    echo ~{default="error" sep="," numbers}
    echo ~{default="success" sep="," missing_numbers}
  >>>

  output {
    String result = read_string(stdout())
    Boolean present = defined(absent)
  }
}
