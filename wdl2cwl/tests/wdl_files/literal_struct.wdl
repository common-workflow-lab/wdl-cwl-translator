version 1.1

struct Foo {
  String one
  String two
}

task literal_struct_test {
  Foo three = { "one": "four", "two": "five" }

  input {
    Foo six = { "one": "seven", "two": "eight"}
  }
  command <<<
    echo ~{three.one} ~{three.two} ~{six.one} ~{six.two}
  >>>

  output {
    String result = read_string(stdout())
  }

}
