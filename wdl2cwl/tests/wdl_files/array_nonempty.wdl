version 1.1
task nonempty_array_test {
  input {
    Array[Int]+ numbers = [ 1, 2, 3]
  }

  command <<<
    echo ~{sep(",", numbers)}
  >>>

  output {
    String result = read_string(stdout())
  }
}
