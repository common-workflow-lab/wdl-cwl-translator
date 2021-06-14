#!/usr/bin/env python3
import sys

from setuptools import setup

#exec(open("cwl_utils/__meta__.py").read())
from setuptools import setup, find_packages

setup(name="wdl2cwl", packages=find_packages())

needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest < 7", "pytest-runner"] if needs_pytest else []

setup(
    install_requires = [
        "cwl-utils",
        "antlr4-python3-runtime",
    ]
)
