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
        ("wdl_files/bcftools_stats.wdl", "oop_cwl_files/bcftools_stats.cwl"),
        ("wdl_files/bcftools_annotate.wdl", "oop_cwl_files/bcftools_annotate.cwl"),
        ("wdl_files/validateOptimus_1.wdl", "oop_cwl_files/validateOptimus_1.cwl"),
        ("wdl_files/validateOptimus_2.wdl", "oop_cwl_files/validateOptimus_2.cwl"),
        ("wdl_files/validateOptimus_3.wdl", "oop_cwl_files/validateOptimus_3.cwl"),
        ("wdl_files/validateOptimus_4.wdl", "oop_cwl_files/validateOptimus_4.cwl"),
        ("wdl_files/validateOptimus_5.wdl", "oop_cwl_files/validateOptimus_5.cwl"),
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
