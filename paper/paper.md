---
title: 'wdl2cwl: A Python tool for converting WDL workflow files to CWL workflows'
tags:
  - Python
  - workflow
  - workflow description language
  - common workflow language
  - Language
date:  August 2022
bibliography: paper.bib
---

# Summary

Computational workflows are used in data analysis, enabling innovation and 
decision-making in various fields of science. There are several competing workflow description languages, 
WDL (Workflow description Language) and CWL (Common Workflow Languaage) are two examples
of popularly used systems. WDL is a way to specify data processiing workflows with 
a human-readable and writeable syntax. WDL makes it straightforward to define 
complex analysis tasks, chain them together in workflows, and parallelize their execution [@wdl].
CWL is a way to describe command line tools and connect 
the together to create workflows. Because CWL is a specification and not a 
specific piece of software, tools and workflows described using CWL are portable acroos
a variety of platforms that support the CWL standard [@cwl:2021].

It can be hard learning and adapting to the various knowledge requirements and
specifications of a workflow language. This complication makes it hard for analysts and scientists 
share their work and collaborate efficiently with others who are use to working with different language.

``wdl2cwl`` is a Python package for converting workflow files written in WDL to workflows in CWL. Python
allows flexibility and the ability to import other Python packages. 
The ``wdl2cwl`` package uses the ``miniwdl`` package [@miniwdl], an open-source, Python package, 
developer toolkit and local runner for running wdl files. The ``wdl2cwl`` package was
designed to provide a class-based and developer-friendly interface for extracting 
the various sections of a wdl file, dynamically create the various objects
needed to create a valid cwl file ( this is done using the ``cwl.utils`` package [@cwl.utils], 
an open-source, Python package for parsing cwl files) and performing error checking 
for valid imports, input declarations, expressions etc.

``wdl2cwl`` was designed to be used by both cwl and wdl professionals and researchers, as well
students in courses on bioinformatics and data analytics. It has already been
used to convert a number of wdl workflows. This tool makes it easier 
to translate your favourite WDL workflows into CWL. This could be because of platform 
dependency issues (like trying to reproduce a WDL workflow on a platform that only allows CWL),
or it could be to eliminate the need for an anylyst who undersnd and writes CWL to learn WDL language
in other to quickly use a WDL workflow or it could be to create a more widely available workflow
that can be used by users of either languages and so much more.

# Acknowledgements

We acknowledge the work by tthe developers of miniwdl, cwl, wdl

# References