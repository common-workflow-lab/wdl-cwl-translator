"""Tests for WDLSourceLine."""

import re

import pytest
import WDL

from ..errors import WDLSourceLine
from ..main import ConversionException
from .util import get_data


def test_wdlsourceline() -> None:
    """Basic test of WDLSourceLine."""
    doc_tree = WDL.load(get_data("wdl_files/minCores.wdl"))
    with pytest.raises(
        ConversionException,
        match=re.escape("minCores.wdl:1:1: SourceLineTest"),
    ):
        raise WDLSourceLine(doc_tree, ConversionException).makeError("SourceLineTest")


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
    doc_tree = WDL.load(get_data("wdl_files/minCores.wdl"))
    with pytest.raises(
        ConversionException,
        match=re.escape("minCores.wdl:3:1: SourceLineTest"),
    ):
        with WDLSourceLine(doc_tree.tasks[0], ConversionException):
            raise WDLSourceLine(
                doc_tree.tasks[0].inputs, ConversionException
            ).makeError("SourceLineTest")
