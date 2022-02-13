cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: read_boolean
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                echo true > true.txt
                echo false > false.txt
                echo True > True.txt
                echo False > False.txt
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        outdirMin: 1024
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: good_true
        type: boolean
        outputBinding:
            loadContents: true
            glob: true.txt
            outputEval: |-
                ${
                var contents = self[0].contents.trim().toLowerCase()
                if (contents == 'true') { return true;}
                if (contents == 'false') { return false;}
                throw "'read_boolean' received neither 'true' nor 'false': " + self[0].contents;
                }
      - id: good_false
        type: boolean
        outputBinding:
            loadContents: true
            glob: false.txt
            outputEval: |-
                ${
                var contents = self[0].contents.trim().toLowerCase()
                if (contents == 'true') { return true;}
                if (contents == 'false') { return false;}
                throw "'read_boolean' received neither 'true' nor 'false': " + self[0].contents;
                }
      - id: mixed_case_true
        type: boolean
        outputBinding:
            loadContents: true
            glob: True.txt
            outputEval: |-
                ${
                var contents = self[0].contents.trim().toLowerCase()
                if (contents == 'true') { return true;}
                if (contents == 'false') { return false;}
                throw "'read_boolean' received neither 'true' nor 'false': " + self[0].contents;
                }
      - id: mixed_case_false
        type: boolean
        outputBinding:
            loadContents: true
            glob: False.txt
            outputEval: |-
                ${
                var contents = self[0].contents.trim().toLowerCase()
                if (contents == 'true') { return true;}
                if (contents == 'false') { return false;}
                throw "'read_boolean' received neither 'true' nor 'false': " + self[0].contents;
                }
  - cwlVersion: v1.2
    id: read_bad_boolean
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                echo 1 > bad-true.txt
                echo 0 > bad-false.txt
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        outdirMin: 1024
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: bad_true
        type: boolean
        outputBinding:
            loadContents: true
            glob: bad-true.txt
            outputEval: |-
                ${
                var contents = self[0].contents.trim().toLowerCase()
                if (contents == 'true') { return true;}
                if (contents == 'false') { return false;}
                throw "'read_boolean' received neither 'true' nor 'false': " + self[0].contents;
                }
      - id: bad_false
        type: boolean
        outputBinding:
            loadContents: true
            glob: bad-false.txt
            outputEval: |-
                ${
                var contents = self[0].contents.trim().toLowerCase()
                if (contents == 'true') { return true;}
                if (contents == 'false') { return false;}
                throw "'read_boolean' received neither 'true' nor 'false': " + self[0].contents;
                }
  - cwlVersion: v1.2
    id: read_dynamic_boolean
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                echo true > $(inputs.filename)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        outdirMin: 1024
    inputs:
      - id: filename
        default: foobar
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: dynamic_true
        type: boolean
        outputBinding:
            loadContents: true
            glob: $(inputs.filename)
            outputEval: |-
                ${
                var contents = self[0].contents.trim().toLowerCase()
                if (contents == 'true') { return true;}
                if (contents == 'false') { return false;}
                throw "'read_boolean' received neither 'true' nor 'false': " + self[0].contents;
                }
