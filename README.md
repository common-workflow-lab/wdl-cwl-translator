# A translator from Workflow Definition Language (WDL) to Common Workflow Language (CWL)

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
    i. https://github.com/broadinstitute/warp
    ii. https://github.com/broadinstitute/viral-pipelines/tree/master/pipes/WDL
    iii. https://github.com/biowdl/tasks
2. Run the workflow using the translator and save the result to a .cwl file at a specified location. ```python wdl2cwl/main.py  path_to_wdl_file -o specified_location```.
3. Add the WDL workflow to ```wdl2cwl/tests/wdl_files``` and the resultant CWL file to ```wdl2cwl/tests/wdl_files```. 
4. Add the paths of the added WDL and CWL files to ```wdl2cwl/tests/test_cwl.py``` as an argument under the @pytest.mark.parametrize() function.
