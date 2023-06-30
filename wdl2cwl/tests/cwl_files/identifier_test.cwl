cwlVersion: v1.2
id: identifer_test
class: CommandLineTool
doc: Javascript reserved words as WDL indentifiers
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |2

            echo $(inputs["do"]) $(inputs["for"]) $(inputs["let"]) $(inputs["new"]) $(inputs["try"]) $(inputs["var"]) $(inputs["case"]) $(inputs["enum"]) $(inputs["eval"])
            echo $(inputs["null"]) $(inputs["void"]) $(inputs["with"]) $(inputs["break"]) $(inputs["catch"]) $(inputs["class"]) $(inputs["const"]) $(inputs["super"])
            echo $(inputs["throw"]) $(inputs["while"]) $(inputs["yield"]) $(inputs["delete"]) $(inputs["export"]) $(inputs["public"]) $(inputs["return"])
            echo $(inputs["static"]) $(inputs["switch"]) $(inputs["typeof"]) $(inputs["default"]) $(inputs["finally"]) $(inputs["package"]) $(inputs["private"])
            echo $(inputs["continue"]) $(inputs["debugger"]) $(inputs["function"]) $(inputs["arguments"]) $(inputs["interface"]) $(inputs["protected"])
            echo $(inputs["implements"]) $(inputs["instanceof"])
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: ResourceRequirement
    outdirMin: 1024
  - class: ToolTimeLimit
    timelimit: 60
inputs:
  - id: do
    type: string
  - id: for
    type: string
  - id: let
    type: string
  - id: new
    type: string
  - id: try
    type: string
  - id: var
    type: string
  - id: case
    type: string
  - id: enum
    type: string
  - id: eval
    type: string
  - id: 'null'
    type: string
  - id: void
    type: string
  - id: with
    type: string
  - id: break
    type: string
  - id: catch
    type: string
  - id: class
    type: string
  - id: const
    type: string
  - id: super
    type: string
  - id: throw
    type: string
  - id: while
    type: string
  - id: yield
    type: string
  - id: delete
    type: string
  - id: export
    type: string
  - id: public
    type: string
  - id: return
    type: string
  - id: static
    type: string
  - id: switch
    type: string
  - id: typeof
    type: string
  - id: default
    type: string
  - id: finally
    type: string
  - id: package
    type: string
  - id: private
    type: string
  - id: continue
    type: string
  - id: debugger
    type: string
  - id: function
    type: string
  - id: arguments
    type: string
  - id: interface
    type: string
  - id: protected
    type: string
  - id: implements
    type: string
  - id: instanceof
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: result
    type: stdout
