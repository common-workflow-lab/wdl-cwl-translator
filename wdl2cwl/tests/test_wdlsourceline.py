"""Tests for WDLSourceLine."""

import re

import pytest
import WDL

from ..errors import WDLSourceLine
from ..main import ConversionException, Converter
from .test_cwl import get_file


def test_wdlsourceline() -> None:
    """Basic test of WDLSourceLine."""
    doc_tree = WDL.load(get_file("wdl_files/minCores.wdl"))
    with pytest.raises(
        ConversionException,
        match=re.escape(
            "minCores.wdl:1:1: Unimplemented type: <class 'WDL.Tree.Document'>"
        ),
    ):
        Converter().load_wdl_objects(doc_tree)  # type: ignore


def test_wdlsourceline_non_wdl() -> None:
    """Using non-miniWDL objects with WDLSourceLine."""
    message = WDLSourceLine({"non-WDL object"}).makeError("Something happened")

    assert message == "Something happened"


def test_wdlsourceline_non_wdl_context_manager() -> None:
    """Test non-miniWDL objects with WDLSourceLine as a context manager."""
    with pytest.raises(ConversionException, match="Something went wrong"):
        with WDLSourceLine({"non-WDL object"}, ConversionException) as test:
            assert test is not None
            raise Exception("Something went wrong")
    assert True


def test_nested_wdlsourceline() -> None:
    """Nested test of WDLSourceLine."""
    doc_tree = WDL.load(get_file("wdl_files/minCores.wdl"))
    with pytest.raises(
        ConversionException,
        match=re.escape(
            "minCores.wdl:1:1: Unimplemented type: <class 'WDL.Tree.Document'>"
        ),
    ):
        with WDLSourceLine(doc_tree.tasks[0], ConversionException):
            Converter().load_wdl_objects(doc_tree)  # type: ignore
