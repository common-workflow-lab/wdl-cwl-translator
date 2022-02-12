# Example based off miniwdl's test_select test function:
# https://github.com/chanzuckerberg/miniwdl/blob/0ee2bddaea4bd2538cd21bd9b7806a6c81068b2f/tests/test_5stdlib.py#L156

version 1.0
task test_select_all {
  input {
    Int one
    Int? two
  }
  command {}
  output {
    Array[Int] first1 = select_all([one, two, 1, select_first([two, one])])
  }

  parameter_meta {
    first1: "amalgamation"
  }
}
