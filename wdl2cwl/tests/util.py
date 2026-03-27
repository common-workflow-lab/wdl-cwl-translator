"""Helper methods for tests."""

import atexit
import os
from contextlib import ExitStack
from importlib.resources import as_file, files
from pathlib import Path


def get_path(filename: str) -> Path:
    """Get the filepath for a given test file."""
    # normalizing path depending on OS or else it will cause problem when joining path
    filename = os.path.normpath(filename)
    file_manager = ExitStack()
    atexit.register(file_manager.close)
    traversable = files("wdl2cwl") / "tests" / filename
    filepath = file_manager.enter_context(as_file(traversable))
    return filepath.resolve()


def get_data(filename: str) -> str:
    return str(get_path(filename))
