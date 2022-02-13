version 1.1
task identifer_test {
  input {
    String do
    String for
    String let
    String new
    String try
    String var
    String case
    String enum
    String eval
    String null
    String void
    String with
    String break
    String catch
    String class
    String const
    String super
    String throw
    String while
    String yield
    String delete
    String export
    String public
    String return
    String static
    String switch
    String typeof
    String default
    String finally
    String package
    String private
    String continue
    String debugger
    String function
    String arguments
    String interface
    String protected
    String implements
    String instanceof
  }

  command <<<
    echo ~{do} ~{for} ~{let} ~{new} ~{try} ~{var} ~{case} ~{enum} ~{eval}
    echo ~{null} ~{void} ~{with} ~{break} ~{catch} ~{class} ~{const} ~{super}
    echo ~{throw} ~{while} ~{yield} ~{delete} ~{export} ~{public} ~{return}
    echo ~{static} ~{switch} ~{typeof} ~{default} ~{finally} ~{package} ~{private}
    echo ~{continue} ~{debugger} ~{function} ~{arguments} ~{interface} ~{protected}
    echo ~{implements} ~{instanceof}
  >>>

  output {
   String result = stdout()
  }

  runtime {
    time_minutes: 1
  }

  meta {
   description: "Javascript reserved words as WDL indentifiers"
  }
}
