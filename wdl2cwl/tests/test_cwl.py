"""Tests for miniwdl."""
import os.path
from tempfile import NamedTemporaryFile
from pathlib import Path
import pytest

from ..main import main


def get_file(path: str) -> str:
    """Get file."""
    return os.path.join(os.path.dirname(__file__), path)


@pytest.mark.parametrize(
    "wdl_path,cwl_path",
    [
        ("wdl_files/bowtie_1.wdl", "cwl_files/bowtie_1.cwl"),
        ("wdl_files/bowtie_2.wdl", "cwl_files/bowtie_2.cwl"),
        ("wdl_files/bowtie_3.wdl", "cwl_files/bowtie_3.cwl"),
        ("wdl_files/bcftools_stats.wdl", "cwl_files/bcftools_stats.cwl"),
        ("wdl_files/bcftools_annotate.wdl", "cwl_files/bcftools_annotate.cwl"),
        ("wdl_files/validateOptimus_1.wdl", "cwl_files/validateOptimus_1.cwl"),
        ("wdl_files/validateOptimus_2.wdl", "cwl_files/validateOptimus_2.cwl"),
        ("wdl_files/validateOptimus_3.wdl", "cwl_files/validateOptimus_3.cwl"),
        ("wdl_files/validateOptimus_4.wdl", "cwl_files/validateOptimus_4.cwl"),
        ("wdl_files/validateOptimus_5.wdl", "cwl_files/validateOptimus_5.cwl"),
        ("wdl_files/deepvariant.wdl", "cwl_files/deepvariant.cwl"),
        ("wdl_files/gatk_1.wdl", "cwl_files/gatk_1.cwl"),
        ("wdl_files/hisat2_1.wdl", "cwl_files/hisat2_1.cwl"),
        ("wdl_files/hisat2_2.wdl", "cwl_files/hisat2_2.cwl"),
        ("wdl_files/hisat2_3.wdl", "cwl_files/hisat2_3.cwl"),
        ("wdl_files/isoseq3.wdl", "cwl_files/isoseq3.cwl"),
        ("wdl_files/minCores.wdl", "cwl_files/minCores.cwl"),
        ("wdl_files/pbmm2.wdl", "cwl_files/pbmm2.cwl"),
        ("wdl_files/rtg_VcfEval.wdl", "cwl_files/rtg_VcfEval.cwl"),
        ("wdl_files/rtg.wdl", "cwl_files/rtg.cwl"),
        ("wdl_files/seqtk_1.wdl", "cwl_files/seqtk_1.cwl"),
        ("wdl_files/seqtk_2.wdl", "cwl_files/seqtk_2.cwl"),
        ("wdl_files/smoove.wdl", "cwl_files/smoove.cwl"),
        ("wdl_files/transcriptclean_1.wdl", "cwl_files/transcriptclean_1.cwl"),
        ("wdl_files/UmiCorrection.wdl", "cwl_files/UmiCorrection.cwl"),
        ("wdl_files/vardict.wdl", "cwl_files/vardict.cwl"),
        ("wdl_files/vt.wdl", "cwl_files/vt.cwl"),
        ("wdl_files/TrimAdapters.wdl", "cwl_files/TrimAdapters.cwl"),
        (
            "wdl_files/read_string_cornercase.wdl",
            "cwl_files/read_string_cornercase.cwl",
        ),
    ],
)
def test_wdls(wdl_path: str, cwl_path: str) -> None:
    """Test WDL to CWL conversion."""
    with open(get_file(cwl_path), encoding="utf-8") as file:
        with NamedTemporaryFile(mode="w", encoding="utf-8") as output:
            main(["-o", output.name, get_file(wdl_path)])
            output.flush()
            assert Path(output.name).read_text(encoding="utf-8") == file.read()


def test_wdl_stdout(capsys) -> None:  # type: ignore
    """Test WDL to CWL conversion using stdout."""
    with open(get_file("cwl_files/bowtie_1.cwl"), encoding="utf-8") as file:
        main([get_file("wdl_files/bowtie_1.wdl")])
        captured = capsys.readouterr()
        assert captured.out == file.read()
        assert captured.err == ""
