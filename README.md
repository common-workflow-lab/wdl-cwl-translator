# A translator from [WDL](https://openwdl.org/) to [CWL v1.2](https://w3id.org/cwl/v1.2/)

[![codecov](https://codecov.io/gh/common-workflow-lab/wdl-cwl-translator/branch/main/graph/badge.svg?token=lvcnJHP1hj)](https://codecov.io/gh/common-workflow-lab/wdl-cwl-translator)

## Background

Workflow Definition Language (WDL) and Common Workflow Language (CWL)
are high-level languages for describing how to run a sequence of
programs to perform a data analysis task.  A workflow consists of a
series of steps that are connected by input/output dependencies.

CWL is the product of community-based open source standards process,
and workflows written in CWL are portable across a number of different
software platforms (e.g. Arvados, Toil, CWL-Airflow, Seven Bridges).  WDL is also
open source, but based largely around a single implementation
(Cromwell), however some workflows that are important to the
bioinformatics community are only maintained in WDL.

The goal of this project is to develop a translator that takes a WDL
workflow and produces an equivalent workflow in CWL.  When executed
with the same input, the translated workflow should produce equivalent
results to the original workflow.  An ideal demonstration of
capability would be to translate the Broad Institute WDL Analysis
Research Pipelines (WARP) Whole Genome Germline Single Sample
workflow, run it on a scale-out, production CWL runner (such as
Arvados or Toil), and show that the results are equivalent.

More background reading on CWL:

- A recent paper: https://arxiv.org/abs/2105.07028 ([Full PDF](https://arxiv.org/pdf/2105.07028.pdf))
- https://www.commonwl.org/user_guide/

This project uses the CWL parser and objects from [`cwl_utils.parser.cwl_v1_2`](https://cwl-utils.readthedocs.io/en/latest/autoapi/cwl_utils/parser/cwl_v1_2/index.html)

[miniwdl](https://github.com/chanzuckerberg/miniwdl) is used for WDL parsing,
and while we target OpenWDL 1.1, earlier versions of (Open)WDL seem to work
thanks to the flexibility of the `miniwdl` parser.

For some discussion comparing the two languages (mainly from the perspective of translating in the other direction, CWL to WDL), see this document:

https://github.com/dnanexus/dxCompiler/blob/main/doc/CWL_v1.2.0_to_WDL_v1.md

## Installation

### Prerequisites

Python 3.7+

### Instructions

These instructions assume a Linux / macOS operating system.

``` shell
git clone https://github.com/common-workflow-lab/wdl-cwl-translator/
cd wdl-cwl-translator
python3 -m venv env
source env/bin/activate
pip install -U pip setuptools wheel
pip install -e .
```

## Usage

``` shell
wdl2cwl path_to_wdl_file
```
To output the CWL version to your terminal/stdout.

``` shell
wdl2cwl path_to_workflow.wdl --output path_to_new_workflow.cwl
```

## Temporary Limitations

WDL features not yet supported
- Advanced [Scatter](https://github.com/common-workflow-lab/wdl-cwl-translator/issues/183)

WDL types not yet supported
- [Map](https://github.com/common-workflow-lab/wdl-cwl-translator/issues/77) types
- [Nested structs](https://github.com/common-workflow-lab/wdl-cwl-translator/issues/158)
- [Pair](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#pairx-y)
- [Object](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-object)

OpenWDL 1.1 standard library functions to be implemented
- [floor](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#int-floorfloat-int-ceilfloat-and-int-roundfloat)
- [min](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-int-minint-int-float-minfloat-float-float-minint-float-float-minfloat-int)
- [max](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-int-maxint-int-float-maxfloat-float-float-maxint-float-float-maxfloat-int)
- [stderr](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#file-stderr)
- [read_map](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#mapstring-string-read_mapstringfile)
- [read_object](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-object-read_objectstringfile)
- [read_objects](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-arrayobject-read_objectsstringfile)
- [read_json](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#r-read_jsonstringfile)
- [write_lines](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#file-write_linesarraystring)
- [write_tsv](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#file-write_tsvarrayarraystring)
- [write_map](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#file-write_mapmapstring-string)
- [write_object](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-file-write_objectobject)
- [write_objects](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-file-write_objectsarrayobject)
- [write_json](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#file-write_jsonx)
- [range](https://github.com/common-workflow-lab/wdl-cwl-translator/issues/176)
- [transpose](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#arrayarrayx-transposearrayarrayx)
- [zip](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#arraypairxy-ziparrayx-arrayy)
- [unzip](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-pairarrayx-arrayy-unziparraypairx-y)
- [cross](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#arraypairxy-crossarrayx-arrayy)
- [prefix](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#arraystring-prefixstring-arrayp)
- [suffix](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-arraystring-suffixstring-arrayp)
- [quote](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-arraystring-quotearrayp)
- [squote](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-arraystring-squotearrayp)
- [as_pairs](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-arraypairp-y-as_pairsmapp-y)
- [as_map](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-mapp-y-as_maparraypairp-y)
- [keys](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-arrayp-keysmapp-y)
- [collect_by_keys](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-mapp-arrayy-collect_by_keyarraypairp-y)

Many of the above are straightforward to implement, but we haven't needed them yet.
So if you are unable to translate a particular WDL document due to lackof a standard library
function, please [open and issue](https://github.com/common-workflow-lab/wdl-cwl-translator/issues/new/choose)
and share your example!

## Incompatibilities possibly requiring manual intervention

1. Dynamic specification of Docker containers.
   As of CWL v1.2, CWL's `DockerRequirement` has no support for dynamic
specifications, only fixed values. If a WDL task has a `runtime.docker` that
references an input with a default value, then `wdl2cwl` does try to copy that
default value to the CWL `DockerRequirement.dockerPull`.

If changing the software container is needed, there are several workarounds:
1. Use workflow runner/engine provided overrides: many CWL runners
(including those based upon the CWL reference runner, `cwltool`) support overriding
requirements at any level at run time. See
https://github.com/common-workflow-language/cwltool#overriding-workflow-requirements-at-load-time
2. Manually override the `DockerRequirement` in `hints` by specifying your own container
at the CWL workflow step level under `requirements`
3. Manually editing the CommandLineTool definition yourself.

### File Localization

(Open)WDL assumes that users will configure localization by placing
input files in the same directory. Descriptions that require this will need
modification before conversion to CWL, as CWL has explicit constructs for
achieving localization (`secondaryFiles`, `InitialWorkDirRequirement`, and/or
explicit staging).

See [this example](https://github.com/mr-c/biowdl_tasks_cwlcompat/commit/0dd4704ec8969e491e6358fe2e8283272cafde21#diff-c76c01f3ca967cdb9c157a75e7fb1a08d0037543b455c2107398601a2f526ebfR45)
for one method using explicit staging of input files in the `command` block to
achieve the localization required by the tool(s) being called.

## Tips for manually improving the CWL outputs

If you are converting a WDL workflow to the CWL format and the original WDL
document is the "source of truth", then one should avoid making manual changes
to the CWL as you will need to maintain those changes as the source WDL
document(s) changes.

Otherwise, for those users looking to convert from WDL to CWL and then continue
to modify the CWL directly, then we have the following advice:

Consider swapping the `wdl2cwl` translation of the WDL tasks for
[community maintained CWL descriptions for popular tools](https://github.com/common-workflow-library/bio-cwl-tools)
when possible. Follow the [instructions on usage](https://github.com/common-workflow-library/bio-cwl-tools#how-to-use-these-descriptions-in-your-own-repository)
and update the `run` line to refer to a local path or a "raw" GitHub URL of the
community-maintained tool description. You may need to adjust a few input names to
match. Of course, we are happy to [receive your enhancements and additional CWL
bio\* tool descriptions](https://github.com/common-workflow-library/bio-cwl-tools#how-to-donate-your-tool-descriptions)!

For the resulting CWL `Workflow` and any CWL `CommandLineTool`s not swapped for
idiomatic CWL descriptions, consider using the following CWL features absent in WDL

### CWL `Workflow` Tips
1. Consider collapsing nested WDL scatters in a single [multi-dimensional CWL scatter](https://www.commonwl.org/v1.2/Workflow.html#Scatter/gather)

### CWL `CommandLineTool` Tips

1. Use `secondaryFiles` instead of implicit file co-localizaton for when you
   have a file and its index(es).
2. Adding [`format`](https://www.commonwl.org/user_guide/16-file-formats/index.html)
   specifiers to input and output `File`s and arrays of `File`s both at the
   `Workflow` and `CommandLineTool` levels. This helps improve the type checking
   of the workflow and anyone wanting to re-use or adapt the individual `CommandLineTool`s.
3. In addition to the `minCores` in ResourceRequirement, consider setting the `maxCores`
   if the tool is known to not benefit from additional cores after a certain amount.
4. Retrieving the actual number of cores allocated via `$(runtime.cores)` to pass
   to your tools.
5. Need the absolute path of the working directory? You can use `$(runtime.outdir)`.
6. Only running a single command and are redirecting a file into in (`my_tool < input_file`)?
   You can change the input to be [`type: stdin`](https://www.commonwl.org/v1.2/CommandLineTool.html#stdin)
   instead of `type: File` and drop the `< input_file` as a shortcut.
7. Specify the underlying tool(s) required beyond a `DockerRequirement` via
   [`SoftwareRequirement`](https://www.commonwl.org/v1.2/CommandLineTool.html#SoftwareRequirement).
   This makes for good documentation, helps give credit to the authors of the tool(s),
   and makes it easier for those who want to run with local software, conda packages,
   and other non-containerized environments.
8. Moving any environment variable settings (`export FOO=bar`) present in the `script.bash` to an
   [`EnvVarRequirement`](https://www.commonwl.org/user_guide/12-env/index.html)
   Be careful, if the `script.bash` runs many commands and the environment variables
   are not set at the beginning, that may be due to them not being appropriate for
   all the commands; so test to confirm that they are safe to move to an `EnvVarRequirement`
   and if you aren't sure, leave them there.
   Per-tool invocations with environment variables like `FOO=bar name_of_tool.pl --option`
   are also a candidate if (1) there are no other tools invoked or (2) they all
   have the same environment variables set or (3) they other tools ignore the
   environment variables.
9. Many WDL `command` sections create output directories and perform other
   "housekeeping" that is not necessary in CWL, like symlinking files to change
   names or otherwise arrange the input files. Output directories that themselves
   don't become a `Directory` type are likely removable. If a
   [specific arrangement of inputs files is needed](https://www.commonwl.org/user_guide/15-staging/index.html),
   or [additional files need to created dynamically](https://www.commonwl.org/user_guide/14-runtime/index.html),
   then consider using `InitialWorkDirRequirement`.
10. Some WDL `command` sections include copying input files to obtain writable versions.
    This can be quite slow on many systems, and from a CWL perspective it is better
    to use [`InitialWorkDirRequirement`](https://www.commonwl.org/user_guide/15-staging/index.html)
    to achieve the same results by marking those inputs as being `writable: true`.
11. Most CWL runners provide methods to monitor task execution in real time,
    so monitoring scripts and other similar techniques can be removed.
12. If the `script.bash` (which comes from the WDL `command` section) meets the
    following criteria, then consider removing it (and the `InitialWorkDirRequirement`
    if otherwise unused) in favor of directly calling your tool using `baseCommand`
    with the name of the executable and any static command line arguments
    and [`arguments`](https://www.commonwl.org/user_guide/08-arguments/index.html)
    with the remaining mix of dynamic and static command line arguments.
    - No `bash` features like `for` loops and `if` statements
    - A single tool is invoked just once; check for un-escaped semicolons `;` which
      means there are multiple commands on a single line.
    - no input/output redirection, use `stdin`, `stdout`, and `stderr` as need be.
    - no pipelining `|`.

## Development

### Running the tests

```
make install-dep
make test  # just the unit tests
make help # to list major makefile targets
make diff_pydocstyle_report # run a diff to show how much changes where made in the docstyle

tox  # all the code checks
tox -l # list of all configured tox environments
tox -e py39-pydocstyle # perform only pydocstyle tests (py39 is the version of the python interpreter you have installed)
```

### Adding Test Cases

1. Find a WDL workflow. Given below are some links that can be used to find these workflows: 
    1. https://github.com/broadinstitute/warp
    2. https://github.com/broadinstitute/viral-pipelines/tree/master/pipes/WDL
    3. https://github.com/biowdl/tasks
2. Use the translator to convert the WDL file and save the result to a .cwl file at a specified location.
   `python wdl2cwl/main.py path_to_wdl_file -o specified_location`.
3. If a problem is encountered during the translation or if the WDL workflow has a feature that has not been
   implemented yet, submit a new issue with a description of the issue at
   https://github.com/common-workflow-lab/wdl-cwl-translator/issues
4. Check whether the outputs of the WDL workflow and the resultant CWL file are equivalent upon giving the required
   inputs. The required inputs can usually be found in the repository of the workflow itself through a keyword search.
   (eg: name of the workflow, name of the tool used in the command section) The WDL workflow can be run using a workflow
   runner like miniwdl. (Refer the documentation https://github.com/chanzuckerberg/miniwdl) The CWL file can be run
   using cwltool (Refer the documentation https://github.com/common-workflow-language/cwltool)
5. Add the WDL workflow to `wdl2cwl/tests/wdl_files` and the resultant CWL file to `wdl2cwl/tests/cwl_files`.
   Include the licence and the original location of the WDL file as a comment at the beginning of the document. 
6. Add the name of the added WDL file to `wdl2cwl/tests/test_cwl.py` as an argument under the
  `@pytest.mark.parametrize()` function.
7. Please run the code checks via `tox`, and fix as many issue as you can on your own. `make format` will fix many things for you!
