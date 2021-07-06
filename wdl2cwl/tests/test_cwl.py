import os.path
import pytest

from .. import main as wdl


def get_file(path: str) -> str:
    return os.path.join(os.path.dirname(__file__), path)


@pytest.mark.parametrize(
    "wdl_path,cwl_path",
    [
        ("wdl_files/smoove.wdl", "cwl_files/smoove.cwl"),
        (
            "wdl_files/CollectQualityYieldMetrics.wdl",
            "cwl_files/CollectQualityYieldMetrics.cwl",
        ),
    ],
)
class TestParametrized:
    """Contains the test functions for WDL to CWL conversion."""
    def test_wdls(self, wdl_path: str, cwl_path: str) -> None:
        """Test WDL to CWL conversion."""
        convertedStr = wdl.convert(get_file(wdl_path))
        testStr = ""
        with open(get_file(cwl_path)) as file:
            testStr = file.read()

        assert convertedStr == testStr

    def test_wdls_patch(
        self, monkeypatch: pytest.MonkeyPatch, wdl_path: str, cwl_path: str
    ) -> None:
        """Test WDL to CWL conversion with patch."""
        monkeypatch.setattr(
            "sys.argv",
            [
                "python",
                get_file(wdl_path),
            ],
        )
        convertedStr = wdl.main()
        testStr = ""
        with open(get_file(cwl_path)) as file:
            testStr = file.read()

        assert convertedStr == testStr
