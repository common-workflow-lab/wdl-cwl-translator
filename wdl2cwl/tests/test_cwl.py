import os.path
import pytest
import pathlib
import filecmp

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
        ("wdl_files/seqtk.wdl", "cwl_files/seqtk.cwl"),
        ("wdl_files/format.wdl", "cwl_files/format.cwl"),
    ],
)
class TestParameterized:
    """Contains the test functions for WDL to CWL conversion."""

    def test_wdls(self, wdl_path: str, cwl_path: str) -> None:
        """Test WDL to CWL conversion."""
        convertedStr = wdl.convert(get_file(wdl_path))
        testStr = ""
        with open(get_file(cwl_path)) as file:
            testStr = file.read()

        assert convertedStr == testStr

    def test_wdls_patch(
        self,
        tmp_path: pathlib.Path,
        monkeypatch: pytest.MonkeyPatch,
        wdl_path: str,
        cwl_path: str,
    ) -> None:
        """Test WDL to CWL conversion with patch."""
        d = tmp_path / "sub"
        d.mkdir()
        p = d / "result.cwl"

        monkeypatch.setattr(
            "sys.argv",
            [
                "python",
                get_file(wdl_path),
                "--output",
                str(p),
            ],
        )

        wdl.main()

        assert filecmp.cmp(p, get_file(cwl_path))
