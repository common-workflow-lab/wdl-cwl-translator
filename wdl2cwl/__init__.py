"""Module to convert (Open)WDL documents to CWL v1.2."""

import logging

_logger = logging.getLogger("wdl2cwl")
_logger.addHandler(logging.StreamHandler())
_logger.setLevel(logging.INFO)
