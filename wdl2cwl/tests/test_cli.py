"""Misc tests for CLI."""

import pytest
import WDL

from ..main import main
from .util import get_data


def test_bad_wdl(capsys: pytest.CaptureFixture[str]) -> None:
    """Confirm miniwdl error reporting for malformed input."""
    with pytest.raises(WDL.Error.SyntaxError):
        main([get_data("parse_error.wdl")])
    captured = capsys.readouterr()
    assert (
        "parse_error.wdl Ln 3 Col 12) Unexpected token Token('$END', '') at line 3, column 12"
        in captured.err
    )
