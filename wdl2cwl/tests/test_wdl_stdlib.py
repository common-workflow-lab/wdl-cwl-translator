"""Tests for the translation of WDL Stdlib functions."""

from tempfile import NamedTemporaryFile

from wdl2cwl.main import convert


def test_select_all() -> None:
    """Tests the ``select_all`` function."""
    # Example based off miniwdl's test_select test function:
    # https://github.com/chanzuckerberg/miniwdl/blob/0ee2bddaea4bd2538cd21bd9b7806a6c81068b2f/tests/test_5stdlib.py#L156
    wdl = R"""
    version 1.0
    task test_select_all {
      input {
        Int one
        Int? two
      }
      command {}
      output {
        Array[Int] first1 = select_all([one, two])
      }
    }
    """
    with NamedTemporaryFile(mode="w", encoding="utf-8") as output:
        output.write(wdl)
        output.flush()
        cwl = convert(output.name)
        assert cwl["outputs"][0] == {
            "id": "first1",
            "type": {"items": "int", "type": "array"},
            "outputBinding": {
                "glob": "$([inputs.one, inputs.two].filter(element => element !== null) )"
            },
        }
