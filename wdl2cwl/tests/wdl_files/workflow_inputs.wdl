version 1.0

struct Foo {
  String one
  String two
}

workflow foo {
  input {
    String first = "one "
    Float third = 3.14159
    Foo fourth = { "one": "fifth", "two": "sixth" }
  }

  call echo {
    input:
      in = fourth.one
  }

  output {
    String first_result = first
    Float third_result = third
    String echo_result = echo.out
  }

  parameter_meta { first: "test coverage example" }
}

task echo {
  input {
    String in
  }
  command <<<
  >>>
  output {
    String out = in
  }
}
