#!/usr/bin/env python3
import sys

from setuptools import find_packages, setup

needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest < 7", "pytest-runner"] if needs_pytest else []

with open("README.md", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="wdl2cwl",
    url="https://github.com/common-workflow-lab/wdl-cwl-translator",
    version="0.0.1",
    author="Dinithi Wickramaratne",
    author_email="diniwick124@gmail.com",
    long_description_content_type="text/markdown",
    packages=find_packages(),
    setup_requires=[] + pytest_runner,
    install_requires=["cwl-utils>=0.32,<0.34", "miniwdl", "regex"],
)
