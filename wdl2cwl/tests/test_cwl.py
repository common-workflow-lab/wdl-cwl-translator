import os.path
import pytest

from .. import main as wdl


def get_file(path: str) -> str:
    return os.path.join(os.path.dirname(__file__), path)

@pytest.mark.parametrize("wdl_path,cwl_path",
    [("wdl_files/smoove.wdl", "cwl_files/smoove.cwl"),
     ("wdl_files/CollectQualityYieldMetrics.wdl", "cwl_files/CollectQualityYieldMetrics.cwl")])
def test_wdls(wdl_path, cwl_path: str) -> None:
    """Test WDL to CWL conversion."""
    convertedStr = wdl.main([get_file(wdl_path)])
    testStr = ""
    with open(get_file(cwl_path)) as file:
        testStr = file.read()

    assert convertedStr == testStr
