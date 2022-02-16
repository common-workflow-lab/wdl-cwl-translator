version 1.0

workflow foo {
  input {
    String first = "one "
    Float third = 3.14159
  }

  output {
    String first_result = first
    Float third_result = third
  }

  parameter_meta { first: "test coverage example" }
}
