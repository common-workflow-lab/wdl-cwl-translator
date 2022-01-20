"""Tests for miniwdl."""
import os.path
from tempfile import NamedTemporaryFile
from pathlib import Path
import pytest

from ..main import main


def get_file(path: str) -> str:
    """Get file."""
    return os.path.join(os.path.dirname(__file__), path)


def test_meta(capsys: pytest.CaptureFixture[str]) -> None:
    """Test meta warning."""
    main([get_file("wdl_files/validateOptimus_3.wdl")])
    assert "----WARNING: SKIPPING META----" in capsys.readouterr().err


# cd wdl2cwl/tests/wdl_files;  for file in *.wdl; do echo "        (\"${file}\"),"; done
@pytest.mark.parametrize(
    "description_name",
    [
        ("bcftools_annotate.wdl"),
        ("bcftools_stats.wdl"),
        ("bowtie_1.wdl"),
        ("bowtie_2.wdl"),
        ("bowtie_3.wdl"),
        ("CollectQualityYieldMetrics.wdl"),
        ("CollectReadgroupBamQualityMetrics.wdl"),
        ("deepvariant.wdl"),
        ("gatk_1.wdl"),
        ("hisat2_1.wdl"),
        ("hisat2_2.wdl"),
        ("hisat2_3.wdl"),
        ("isoseq3.wdl"),
        # ("length_function_error.wdl"), produces an error on purpose
        ("minCores.wdl"),
        ("pbmm2.wdl"),
        ("read_string_cornercase.wdl"),
        ("rtg_VcfEval.wdl"),
        ("rtg.wdl"),
        ("seqtk_1.wdl"),
        ("seqtk_2.wdl"),
        ("smoove.wdl"),
        ("transcriptclean_1.wdl"),
        ("TrimAdapters.wdl"),
        ("UmiCorrection.wdl"),
        ("validateOptimus_1.wdl"),
        ("validateOptimus_2.wdl"),
        ("validateOptimus_3.wdl"),
        ("validateOptimus_4.wdl"),
        ("validateOptimus_5.wdl"),
        ("vardict.wdl"),
        ("vt.wdl"),
    ],
)
def test_wdls(description_name: str) -> None:
    """Test WDL to CWL conversion."""
    with open(
        get_file(f"cwl_files/{description_name[:-3]}cwl"), encoding="utf-8"
    ) as file:
        with NamedTemporaryFile(mode="w", encoding="utf-8") as output:
            main(["-o", output.name, get_file(f"wdl_files/{description_name}")])
            output.flush()
            assert Path(output.name).read_text(encoding="utf-8") == file.read()


def test_wdl_stdout(capsys) -> None:  # type: ignore
    """Test WDL to CWL conversion using stdout."""
    with open(get_file("cwl_files/bowtie_1.cwl"), encoding="utf-8") as file:
        main([get_file("wdl_files/bowtie_1.wdl")])
        captured = capsys.readouterr()
        assert captured.out == file.read()
        assert captured.err == (
            "----WARNING: SKIPPING REQUIREMENT memory----\n"
            "----WARNING: SKIPPING REQUIREMENT disks----\n"
            "----WARNING: SKIPPING REQUIREMENT time_minutes----\n"
            "----WARNING: SKIPPING META----\n"
        )
