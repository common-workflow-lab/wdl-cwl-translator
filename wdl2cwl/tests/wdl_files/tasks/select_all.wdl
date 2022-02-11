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
}