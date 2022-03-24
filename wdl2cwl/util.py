"""Temporary module for WDL functions that do not belong to the converter."""

from dataclasses import dataclass, field
from typing import List, Set, Tuple

import regex  # type: ignore

from wdl2cwl.errors import ConversionException

# is a combination of https://github.com/tc39/proposal-regexp-unicode-property-escapes#other-examples
# and regex at the bottom of https://stackoverflow.com/a/9392578
# double checked against https://262.ecma-international.org/5.1/#sec-7.6
# eval is not on the official list of reserved words, but it is a built-in function
valid_js_identifier = regex.compile(
    r"^(?!(?:do|if|in|for|let|new|try|var|case|else|enum|eval|null|this|true|"
    r"void|with|break|catch|class|const|false|super|throw|while|yield|delete|export|"
    r"import|public|return|static|switch|typeof|default|extends|finally|package|"
    r"private|continue|debugger|function|arguments|interface|protected|implements|"
    r"instanceof)$)(?:[$_\p{ID_Start}])(?:[$_\u200C\u200D\p{ID_Continue}])*$"
)


def get_input(input_name: str) -> str:
    """Produce a concise, valid CWL expr/param reference lookup string for a given input name."""
    if valid_js_identifier.match(input_name):
        return f"inputs.{input_name}"
    return f'inputs["{input_name}"]'


def get_mem_in_bytes(unit: str) -> Tuple[int, int]:
    """
    Determine the value of a memory unit in bytes.

    Returns the base and exponent, ready for stringifying or evaluation
    """
    """Determine the value of a memory unit in bytes."""
    if unit == "KiB" or unit == "Ki":
        return 1024, 1
    elif unit == "MiB" or unit == "Mi":
        return 1024, 2
    elif unit == "GiB" or unit == "Gi":
        return 1024, 3
    elif unit == "TiB" or unit == "Ti":
        return 1024, 4
    elif unit == "B":
        return 1024, 0
    elif unit == "KB" or unit == "K":
        return 1000, 1
    elif unit == "MB" or unit == "M":
        return 1000, 2
    elif unit == "GB" or unit == "G":
        return 1000, 3
    elif unit == "TB" or unit == "T":
        return 1000, 4
    raise ConversionException(f"Invalid memory unit: ${unit}")


def nice_quote(value: str) -> str:
    """Surround string with quotes, with minimal escaping."""
    single = "'" in value
    double = '"' in value
    if not double:
        return f'"{value}"'
    if not single and double:
        return f"'{value}'"
    # single and double quotes found
    return '"' + value.replace('"', r"\"") + '"'


@dataclass
class ConversionContext:
    """Converstion context, including variables used throughout the conversion process."""

    non_static_values: Set[str] = field(default_factory=set)
    optional_cwl_null: Set[str] = field(default_factory=set)
    scatter_names: List[str] = field(default_factory=list)


__all__ = ["get_input", "get_mem_in_bytes", "nice_quote", "ConversionContext"]
