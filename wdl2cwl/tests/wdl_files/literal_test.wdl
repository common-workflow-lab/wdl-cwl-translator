version 1.1

task literal {
  input {
   Boolean first = true
   Int second = 42
   Float third = 6.022
   String fourth = "hoopla"
   File fifth = "../../../README.md"
   Array[Boolean] sixth = [true, false]
   Array[Int] seventh = [42, 23]
   Array[Float] eighth = [6.022, 10.0, 23.0]
   Array[String] nineth = ["Hello", "World"]
   Array[File] tenth = ["../../../README.md", "../../../LICENSE"]
  }

  command <<<
   echo ~{sep(" ", tenth)}
  >>>

  output {
   String result = "~{first} ~{second} ~{third} ~{fourth} ~{basename(fifth)} ~{sep(',', sixth)} ~{sep(',', seventh)} ~{sep(',', eighth)} ~{sep(',', nineth)} ~{basename(tenth[0])},~{basename(tenth[1])}"
  }
}

