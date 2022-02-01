"""Tests for miniwdl."""
import os.path
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, NamedTuple, Union, cast

import pytest

from ..main import ConversionException, Converter, main


def get_file(path: str) -> str:
    """Get file."""
    return os.path.join(os.path.dirname(__file__), path)


def test_meta(caplog: pytest.LogCaptureFixture) -> None:
    """Test meta warning."""
    main([get_file("wdl_files/TrimAdapters.wdl")])
    assert "Skipping meta" in caplog.text


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
        ("rtg_2_tasks.wdl"),
        ("bcftools.wdl"),
        ("align_and_count.wdl"),
        ("BuildCembaReferences.wdl"),
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


def test_wdl_stdout(
    caplog: pytest.LogCaptureFixture, capsys: pytest.CaptureFixture[str]
) -> None:
    """Test WDL to CWL conversion using stdout."""
    with open(get_file("cwl_files/bowtie_1.cwl"), encoding="utf-8") as file:
        main([get_file("wdl_files/bowtie_1.wdl")])
        captured = capsys.readouterr()
        log = caplog.text
        assert captured.out == file.read()
        assert "Skiping requirement: memory" in log
        assert "Skiping requirement: disks" in log
        assert "Skiping requirement: time_minutes" in log


class TestObject(NamedTuple):
    """Test object for creating WDL.Expr.String."""

    value: Union[Any, None]
    literal: Union[Any, None]


testdata = [
    ("20 MB", 19.073486328125),
    ("400 G", 381469.7265625),
    ("350 K", 0.3337860107421875),
    ("5 T", 4768371.58203125),
    ("7000 B", 0.00667572021484375),
    ("90000 KiB", 87.890625),
    ("15000 GiB", 15360000),
    ("5 TiB", 5242880),
    ("6000 MiB", 6000),
    ("1 Altuve", ConversionException),
]


@pytest.mark.parametrize("memory_runtime, expected_memory_or_error", testdata)
def test_get_memory_literal(
    memory_runtime: str, expected_memory_or_error: Union[float, int, Exception]
) -> None:
    """Test get_memory_literal conditional statements."""
    b = TestObject(value=memory_runtime, literal=None)
    a = TestObject(value=None, literal=b)
    if isinstance(expected_memory_or_error, (float, int)):
        assert Converter().get_memory_literal(a) == expected_memory_or_error  # type: ignore
    else:
        with pytest.raises(expected_memory_or_error):  # type: ignore
            Converter().get_memory_literal(a)  # type: ignore
