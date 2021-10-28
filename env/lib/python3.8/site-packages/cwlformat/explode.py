# Copyright (c) 2020 Seven Bridges
# Given a cwl dict, if it is a workflow, split out any inlined steps into their own files

from typing import List
import pathlib
import sys

import ruamel.yaml

from .version import __version__
from .formatter import stringify_dict, leading_comment_lines


yaml = ruamel.yaml.YAML()


class CWLProcess:
    def __init__(self, cwl: dict, file_path: pathlib.Path):
        self.cwl = cwl
        self.file_path = file_path

    def __str__(self):
        return stringify_dict(self.cwl)

    def save(self):
        self.file_path.parent.mkdir(parents=True, exist_ok=True)
        self.file_path.write_text(leading_comment_lines("") + stringify_dict(self.cwl))


def explode(cwl: CWLProcess) -> List[CWLProcess]:
    _processes = [cwl]
    _cwl = cwl.cwl
    if _cwl.get("class") == "Workflow":
        _cwl_steps = _cwl.get("steps", {})
        _is_dict = isinstance(_cwl_steps, dict)
        for _k, _step in (_cwl_steps.items() if _is_dict else enumerate(_cwl_steps)):
            _step_id = _k if _is_dict else _step.get("id")
            if _step_id is not None:
                _run = _step.get("run")
                if isinstance(_run, dict):
                    step_path = \
                        cwl.file_path.parent / \
                        (cwl.file_path.name + ".steps") / \
                        (_step_id + ".cwl")
                    _step["run"] = str(step_path.relative_to(cwl.file_path.parent))
                    _processes += explode(CWLProcess(_run, step_path))

    return _processes


def main():
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"Rabix/cwl-explode v{__version__}\n"
                    "Explodes CWL workflow with inlined steps")
    parser.add_argument("cwlfile")
    parser.add_argument("outname")
    args = parser.parse_args()

    fp = pathlib.Path(args.cwlfile).absolute()
    as_dict = yaml.load(fp.read_text())
    fp_out = pathlib.Path(args.outname).absolute()
    for n, exploded in enumerate(explode(CWLProcess(as_dict, fp_out))):
        sys.stderr.write(f"{n + 1}: {exploded.file_path.relative_to(fp.parent)}\n")
        exploded.save()


if __name__ == "__main__":
    main()
