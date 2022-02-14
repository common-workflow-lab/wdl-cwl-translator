"""ContextManager to help with error reporting in miniWDL parsed files."""

import re
from types import TracebackType
from typing import Any, Callable, Optional, Type, cast

import WDL
from WDL import SourcePosition

# Inspired by https://github.com/common-workflow-language/schema_salad/blob/661fb0fa8c745ed70253dda93bd12002007f6b33/schema_salad/sourceline.py#L232


lineno_re = re.compile("^(.*?:[0-9]+:[0-9]+: )(( *)(.*))")


class WDLSourceLine:
    """Contextmanager wrap exceptions with WDL source file locations."""

    def __init__(
        self,
        item: Any,
        raise_type: Callable[[str], Any] = str,
    ):
        """Which item and exception type to raise."""
        self.item = item
        self.raise_type = raise_type

    def __enter__(self) -> "WDLSourceLine":
        """Enter the context."""
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_value: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> None:
        """
        Process the exit from the context.

        If there was an exception, wrap it in the raise_type.
        """
        if not exc_value:
            return
        raise self.makeError(str(exc_value)) from exc_value

    def makeLead(self) -> str:
        """Calculate the error message prefix."""
        pos: SourcePosition = cast(SourcePosition, self.item.pos)
        return f"{pos.uri}:{pos.line}:{pos.column}:"

    def makeError(self, msg: str) -> Any:
        """Add the source info to the msg and instantiate the raise_type with it."""
        if not isinstance(
            self.item,
            (
                WDL.Error.SourceNode,
                WDL.Tree.SourceComment,
                WDL.Tree.DocImport,
                WDL.Type.Base,
            ),
        ):
            return self.raise_type(msg)
        errs = []
        lead = self.makeLead()
        for m in msg.splitlines():
            if bool(lineno_re.match(m)):
                errs.append(m)
            else:
                errs.append(f"{lead} {m}")
        return self.raise_type("\n".join(errs))
