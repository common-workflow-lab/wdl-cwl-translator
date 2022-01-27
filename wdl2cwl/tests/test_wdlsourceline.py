"""Tests for WDLSourceLine."""

import re

import pytest
import WDL

from ..main import ConversionException, Converter
from .test_cwl import get_file


def test_wdlsourceline() -> None:
    """Basic test of WDLSourceLine."""
    doc_tree = WDL.load(get_file("wdl_files/minCores.wdl"))
    with pytest.raises(
        ConversionException,
        match=re.escape(
            "minCores.wdl:1:1:Unimplemented type: <class 'WDL.Tree.Document'>"
        ),
    ):
        Converter().load_wdl_objects(doc_tree)  # type: ignore
