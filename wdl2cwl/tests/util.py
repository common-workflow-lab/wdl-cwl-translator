"""Helper methods for tests."""

import os

from pkg_resources import Requirement, ResolutionError, resource_filename


def get_data(filename: str) -> str:
    """Safely get a data file even if we are installed or running from an archive."""
    filename = os.path.normpath(filename)  # normalizing path depending on OS
    # or else it will cause problem when joining path
    filepath = None
    try:
        filepath = resource_filename(Requirement.parse("wdl2cwl"), filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), filename)
    return filepath
