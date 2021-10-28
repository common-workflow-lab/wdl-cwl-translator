# Copyright (c) 2020 Seven Bridges

from typing import Union
import sys
import pathlib
try:
    from importlib.resources import read_text
except ImportError:
    # Python 3.6 fallback: https://importlib-resources.readthedocs.io/en/latest/
    from importlib_resources import read_text

import ruamel.yaml
from ruamel.yaml import scalarstring
from ruamel.yaml.compat import StringIO
from ruamel.yaml.comments import CommentedMap

from cwlformat.version import __version__

yaml = ruamel.yaml.YAML()
yaml.indent(mapping=2, sequence=2, offset=0)
Literal = ruamel.yaml.scalarstring.LiteralScalarString

key_order_dict = yaml.load(read_text("cwlformat", "keyorder.yml"))

hash_bang = "#!/usr/bin/env cwl-runner\n\n"
hash_bang_pre = "#!/usr/bin/env "


def leading_comment_lines(raw_cwl: str):

    top_comment = []

    if len(raw_cwl) > 0 and raw_cwl.lstrip()[0] != "{":
        for _line in raw_cwl.splitlines(keepends=True):
            line = _line.strip()
            if line == "" or line[0] == "#":
                top_comment += [_line]
            else:
                break

    if len(top_comment) == 0 or not top_comment[0].startswith(hash_bang_pre):
        top_comment = [hash_bang] + top_comment

    return "".join(top_comment)


def format_node(cwl: Union[dict, list, str], node_path=None):
    if isinstance(cwl, str):
        if len(cwl) > 80:
            return Literal(cwl)
        else:
            return cwl

    elif isinstance(cwl, dict):
        _fmt_cwl = CommentedMap([
            (k, format_node(v, node_path + [k])) for k, v in reorder_node(cwl, node_path)])
        if _fmt_cwl.get("class") in ["CommandLineTool", "ExpressionTool", "Workflow"]:
            add_space_between_main_sections(_fmt_cwl)
        return _fmt_cwl

    elif isinstance(cwl, list):
        return [format_node(v, node_path) for v in cwl]

    else:
        return cwl


def add_space_between_main_sections(cwl: CommentedMap):
    for k in cwl.keys():
        if k in ["inputs", "outputs", "steps", "requirements", "hints", "baseCommand"]:
            cwl.yaml_set_comment_before_after_key(key=k, before="\n")


def reorder_node(cwl: dict, node_path: list) -> dict:
    known_key_order = key_order_dict.get(
        infer_type(cwl, node_path), key_order_dict["generic-ordering"])
    extra_keys = sorted(set(cwl.keys()) - set(known_key_order))

    for k in known_key_order + extra_keys:
        if k in cwl:
            yield k, cwl[k]


def infer_type(cwl: dict, node_path: list):
    if "class" in cwl:
        return cwl["class"]
    else:
        return "generic-ordering"


def cwl_format(raw_cwl: str) -> str:
    as_dict = yaml.load(raw_cwl)
    return leading_comment_lines(raw_cwl) + stringify_dict(as_dict)


def stringify_dict(as_dict: dict) -> str:
    as_dict = format_node(as_dict, node_path=[])
    stream = StringIO()
    yaml.dump(as_dict, stream)
    return stream.getvalue()


def main():
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"Rabix/cwl-format v{__version__}\n"
                    "A very opinionated code formatter for CWL")
    parser.add_argument("cwlfile")
    parser.add_argument("--inplace", action="store_true",
                        help="Instead of writing formatted code to stdout, overwrite original file")
    args = parser.parse_args()
    fp = pathlib.Path(args.cwlfile)
    formatted = cwl_format(fp.read_text())
    if args.inplace:
        fp.write_text(formatted)
    else:
        sys.stdout.write(formatted)


if __name__ == "__main__":
    main()
