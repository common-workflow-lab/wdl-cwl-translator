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

# For Outreachy contributors

Applicants can demonstrate their technical ability to contribute to
the project by doing the following:

1. Get the ANTLR grammar for WDL 1.1 from https://github.com/openwdl/wdl
2. Generate Python3 parser code from the grammar using the ANTLR tool
https://www.antlr.org/index.html The command for this will be something like:
`antlr4 -Dlanguage=Python -o . grammar.g4`
3. Using the generated parser code, write a simple program that prints
out the names of the ‘tasks’ from
https://github.com/broadinstitute/warp/blob/develop/tasks/broad/Utilities.wdl
4. Using
   [cwl-utils](https://github.com/common-workflow-language/cwl-utils)
   write a simple program that starts with `cwl_utils.parser_v1_2.CommandLineTool` and constructs a CWL data structure that corresponds to this file:
   https://github.com/common-workflow-language/cwl-utils/blob/main/testdata/md5sum.cwl
   then saves it to a `.cwl` file.
5. Make a pull request against this repository with your sample programs
