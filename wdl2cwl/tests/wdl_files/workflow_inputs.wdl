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
      in = fourth.one + ".suffix",
      other = fourth.two,
      echo = true
  }

  output {
    String first_result = first
    Float third_result = third
    String echo_out = echo.out
    String echo_result = echo.result
  }

  parameter_meta { first: "test coverage example" }
}

task echo {
  input {
    String in
    String other
    Boolean echo
  }
  command <<<
    ~{if(echo) then "echo " + in + " " + other else ""}
  >>>
  output {
    String out = in
    String result = read_string(stdout())
  }
}
