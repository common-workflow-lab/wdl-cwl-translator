# A translator from [OpenWDL v1.1](https://github.com/openwdl/wdl/tree/main/versions/1.1) to [CWL v1.2](https://w3id.org/cwl/v1.2/)

Workflow Definition Language (WDL) and Common Workflow Language (CWL)
are high-level languages for describing how to run a sequence of
programs to perform a data analysis task.  A workflow consists of a
series of steps that are connected by input/output dependencies.

CWL is the product of community-based open source standards process,
and workflows written in CWL are portable across a number of different
software platforms (e.g. Arvados, Toil, CWL-Airflow).  WDL is also
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

There is a previous proof of concept called
[wdl2cwl](https://github.com/common-workflow-lab/wdl2cwl), based on
much older versions of WDL and CWL.  It should be used for
inspiration, however we expect this project will involve writing a new
program that uses the ANTLR grammar for WDL.

For some discussion comparing the two languages (mainly from the perspective of translating in the other direction, CWL to WDL), see this document:

https://github.com/dnanexus/dxCompiler/blob/main/doc/CWL_v1.2.0_to_WDL_v1.md

### Adding Test Cases

1. Find a WDL workflow. Given below are some links that can be used to find these workflows: 
    1. https://github.com/broadinstitute/warp
    2. https://github.com/broadinstitute/viral-pipelines/tree/master/pipes/WDL
    3. https://github.com/biowdl/tasks
2. Use the translator to convert the WDL file and save the result to a .cwl file at a specified location. ```python wdl2cwl/main.py  path_to_wdl_file -o specified_location```.
3. If a problem is encountered during the translation or if the WDL workflow has a feature that has not been implemented yet, submit a new issue with a description of the issue at https://github.com/common-workflow-lab/wdl-cwl-translator/issues
4. Check whether the outputs of the WDL workflow and the resultant CWL file are equivalent upon giving the required inputs. The required inputs can usually be found in the repository of the workflow itself through a keyword search. (eg: name of the workflow, name of the tool used in the command section) The WDL workflow can be run using a workflow runner like miniwdl. (Refer the documentation https://github.com/chanzuckerberg/miniwdl) The CWL file can be run using cwltool (Refer the documentation https://github.com/common-workflow-language/cwltool)
5. Add the WDL workflow to ```wdl2cwl/tests/wdl_files``` and the resultant CWL file to ```wdl2cwl/tests/wdl_files```. Include the licence and the original location of the WDL file as a comment at the beginning of the document. 
6. Add the paths of the added WDL and CWL files to ```wdl2cwl/tests/test_cwl.py``` as an argument under the @pytest.mark.parametrize() function.
