"""Tests for miniwdl."""
import os.path
import pytest
import pathlib
import filecmp

from .. import main_oop as wdl


def get_file(path: str) -> str:
    """Get file."""
    return os.path.join(os.path.dirname(__file__), path)


@pytest.mark.parametrize(
    "wdl_path,cwl_path",
    [
        ("wdl_files/bowtie_1.wdl", "oop_cwl_files/bowtie_1.cwl"),
    ],
)
class TestParameterized:
    """Contains the test functions for WDL to CWL conversion."""

    def test_wdls(self, wdl_path: str, cwl_path: str) -> None:
        """Test WDL to CWL conversion."""
        convertedStr = wdl.Converter.load_wdl_tree(get_file(wdl_path))
        testStr = ""
        with open(get_file(cwl_path)) as file:
            testStr = file.read()

        assert convertedStr == testStr
