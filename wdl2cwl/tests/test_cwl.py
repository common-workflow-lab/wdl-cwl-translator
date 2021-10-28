import os.path
import pytest
import pathlib
import filecmp

from .. import main as wdl


def get_file(path: str) -> str:
    return os.path.join(os.path.dirname(__file__), path)


def test_meta(capsys: pytest.CaptureFixture[str]) -> None:
    """Test meta warning."""
    wdl.convert(get_file("wdl_files/UmiCorrection.wdl"))
    assert "----WARNING: SKIPPING META----" in capsys.readouterr().err


@pytest.mark.parametrize(
    "wdl_path,cwl_path",
    [
        ("wdl_files/smoove.wdl", "cwl_files/smoove.cwl"),
        (
            "wdl_files/CollectQualityYieldMetrics.wdl",
            "cwl_files/CollectQualityYieldMetrics.cwl",
        ),
        ("wdl_files/seqtk_1.wdl", "cwl_files/seqtk_1.cwl"),
        ("wdl_files/seqtk_2.wdl", "cwl_files/seqtk_2.cwl"),
        ("wdl_files/hisat2_1.wdl", "cwl_files/hisat2_1.cwl"),
        ("wdl_files/hisat2_2.wdl", "cwl_files/hisat2_2.cwl"),
        (
            "wdl_files/CollectReadgroupBamQualityMetrics.wdl",
            "cwl_files/CollectReadgroupBamQualityMetrics.cwl",
        ),
        ("wdl_files/rtg.wdl", "cwl_files/rtg.cwl"),
        ("wdl_files/bowtie_1.wdl", "cwl_files/bowtie_1.cwl"),
        ("wdl_files/bowtie_2.wdl", "cwl_files/bowtie_2.cwl"),
        ("wdl_files/bowtie_3.wdl", "cwl_files/bowtie_3.cwl"),
        ("wdl_files/UmiCorrection.wdl", "cwl_files/UmiCorrection.cwl"),
        ("wdl_files/vardict.wdl", "cwl_files/vardict.cwl"),
        ("wdl_files/rtg.wdl", "cwl_files/rtg.cwl"),
        ("wdl_files/rtg_2.wdl", "cwl_files/rtg_2.cwl"),
        ("wdl_files/pbmm2.wdl", "cwl_files/pbmm2.cwl"),
        ("wdl_files/isoseq3.wdl", "cwl_files/isoseq3.cwl"),
        ("wdl_files/TrimAdapters.wdl", "cwl_files/TrimAdapters.cwl"),
        ("wdl_files/vt.wdl", "cwl_files/vt.cwl"),
        ("wdl_files/transcriptclean_1.wdl", "cwl_files/transcriptclean_1.cwl"),
        ("wdl_files/gatk_1.wdl", "cwl_files/gatk_1.cwl"),
        ("wdl_files/deepvariant.wdl", "cwl_files/deepvariant.cwl"),
        (
            "wdl_files/validateOptimus_1.wdl",
            "cwl_files/validateOptimus_1.cwl",
        ),
        (
            "wdl_files/validateOptimus_2.wdl",
            "cwl_files/validateOptimus_2.cwl",
        ),
        (
            "wdl_files/validateOptimus_3.wdl",
            "cwl_files/validateOptimus_3.cwl",
        ),
        (
            "wdl_files/validateOptimus_4.wdl",
            "cwl_files/validateOptimus_4.cwl",
        ),
        (
            "wdl_files/validateOptimus_5.wdl",
            "cwl_files/validateOptimus_5.cwl",
        ),
        (
            "wdl_files/bcftools_annotate.wdl",
            "cwl_files/bcftools_annotate.cwl",
        ),
        (
            "wdl_files/bcftools_stats.wdl",
            "cwl_files/bcftools_stats.cwl",
        ),
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


@pytest.mark.parametrize(
    "length_func_error_wdl_file, expected_error_message",
    [
        (
            "wdl_files/length_function_error.wdl",
            "Length function without the if...else keywords is currently not supported",
        ),
    ],
)
def test_length_function_error(
    length_func_error_wdl_file: str, expected_error_message: str
) -> None:
    """Test error message raised when length function does not include the if...else keywords."""
    with pytest.raises(ValueError) as info:
        test_run = wdl.convert(get_file(length_func_error_wdl_file))
    assert str(info.value) == expected_error_message
