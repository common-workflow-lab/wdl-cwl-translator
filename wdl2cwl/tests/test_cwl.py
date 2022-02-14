"""Tests for miniwdl."""
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, NamedTuple, Union

import pytest

from ..main import ConversionException, Converter, main
from .util import get_data


def test_meta(caplog: pytest.LogCaptureFixture) -> None:
    """Test meta warning."""
    main([get_data("wdl_files/InternalTasks.wdl")])
    assert "Skipping meta" in caplog.text


# cd wdl2cwl/tests/wdl_files;  for file in *.wdl; do echo "        (\"${file}\"),"; done
@pytest.mark.parametrize(
    "description_name",
    [
        ("bcftools.wdl"),
        ("bowtie.wdl"),
        ("Qc.wdl"),
        ("deepvariant.wdl"),
        ("gatk.wdl"),
        ("hisat2.wdl"),
        ("isoseq3.wdl"),
        ("minCores.wdl"),
        ("pbmm2.wdl"),
        ("read_string_cornercase.wdl"),
        ("read_boolean.wdl"),
        ("rtg.wdl"),
        ("seqtk.wdl"),
        ("smoove.wdl"),
        ("transcriptclean.wdl"),
        ("TrimAdapters.wdl"),
        ("UmiCorrection.wdl"),
        ("validateOptimus_1.wdl"),
        ("validateOptimus_2.wdl"),
        ("validateOptimus_3.wdl"),
        ("validateOptimus_4.wdl"),
        ("validateOptimus_5.wdl"),
        ("vardict.wdl"),
        ("vt.wdl"),
        ("align_and_count.wdl"),
        ("BuildCembaReferences.wdl"),
        ("bwa.wdl"),
        ("InternalTasks.wdl"),
        ("flatten.wdl"),
        ("select_all_1.wdl"),
        ("merge_svs.wdl"),
        ("literal_test.wdl"),
        ("align_and_count_multiple_report.wdl"),
        ("identifier_test.wdl"),

    ],
)
def test_wdls(description_name: str) -> None:
    """Test WDL to CWL conversion."""
    with open(
        get_data(f"cwl_files/{description_name[:-3]}cwl"), encoding="utf-8"
    ) as file:
        with NamedTemporaryFile(mode="w", encoding="utf-8") as output:
            main(["-o", output.name, get_data(f"wdl_files/{description_name}")])
            output.flush()
            assert Path(output.name).read_text(encoding="utf-8") == file.read()


def test_wdl_stdout(
    caplog: pytest.LogCaptureFixture, capsys: pytest.CaptureFixture[str]
) -> None:
    """Test WDL to CWL conversion using stdout."""
    with open(get_data("cwl_files/bowtie.cwl"), encoding="utf-8") as file:
        main([get_data("wdl_files/bowtie.wdl")])
        captured = capsys.readouterr()
        log = caplog.text
        assert captured.out == file.read()


def test_wdl_url(
    caplog: pytest.LogCaptureFixture, capsys: pytest.CaptureFixture[str]
) -> None:
    """Test WDL to CWL conversion using a HTTPS URL."""
    with open(get_data("cwl_files/bowtie.cwl"), encoding="utf-8") as file:
        main(
            [
                "https://github.com/biowdl/tasks/raw/c5cfd2f5acc2ff729987b86d38b29af046677fdc/bowtie.wdl"
            ]
        )
        captured = capsys.readouterr()
        log = caplog.text
        assert captured.out == file.read()


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
