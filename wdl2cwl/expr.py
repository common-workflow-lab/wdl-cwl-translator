"""WDL Expressions (literal values, arithmetic, comparison, conditional, string interpolation, array, map, and functions)."""

import re
from typing import Any, cast, Optional, Union, Tuple

import WDL
from wdl2cwl.errors import WDLSourceLine, ConversionException
from wdl2cwl.util import get_input, nice_quote, ConversionContext


def get_expr_ifthenelse(
    wdl_ifthenelse: WDL.Expr.IfThenElse, ctx: ConversionContext
) -> str:
    """Translate WDL IfThenElse Expressions."""
    condition = get_expr(wdl_ifthenelse.condition, ctx)
    if_true = get_expr(wdl_ifthenelse.consequent, ctx)
    if_false = get_expr(wdl_ifthenelse.alternative, ctx)
    return f"{condition} ? {if_true} : {if_false}"


def translate_wdl_placeholder(
    wdl_placeholder: WDL.Expr.Placeholder, ctx: ConversionContext
) -> str:
    """Translate WDL Expr Placeholder to a valid CWL expression."""
    expr = wdl_placeholder.expr
    placeholder_expr = get_expr(expr, ctx)
    options = wdl_placeholder.options
    if options:
        if "true" in options:
            true_str = nice_quote(options["true"])
            false_str = nice_quote(options["false"])
            test_str = f"{placeholder_expr} ? {true_str} : {false_str}"
            is_optional = False
            if isinstance(expr, WDL.Expr.Get):
                is_optional = expr.type.optional
            elif isinstance(expr, WDL.Expr.Apply):
                is_optional = (
                    expr.arguments[0].type.optional
                    and expr.function_name != "defined"  # optimization
                )
            if not is_optional:
                return test_str
            else:
                if "default" in options:
                    return (
                        f"{placeholder_expr} === null ? "
                        f"{nice_quote(options['default'])} : {test_str}"
                    )
                return f'{placeholder_expr} === null ? "" : {test_str}'
        elif "sep" in options:
            separator = options["sep"]
            assert isinstance(expr.type, WDL.Type.Array)
            item_type = expr.type.item_type
            if isinstance(item_type, WDL.Type.File):
                pl_holder_str = (
                    placeholder_expr + ".map(function(el) {return el.path})"
                    f'.join("{separator}")'
                )
            else:
                pl_holder_str = f'{placeholder_expr}.join("{separator}")'
            if "default" in options and (expr.type.optional or item_type.optional):
                return (
                    f"{placeholder_expr} === null ? "
                    f"{nice_quote(options['default'])} : {pl_holder_str}"
                )
            return pl_holder_str
        # options must contain only "default", no "sep" or "true"/"false"
        return (
            f"{placeholder_expr} === null ? "
            f"{nice_quote(options['default'])} : {placeholder_expr}"
        )
    return placeholder_expr


def get_expr_string(wdl_expr_string: WDL.Expr.String, ctx: ConversionContext) -> str:
    """Translate WDL String Expressions."""
    if wdl_expr_string.literal is not None:
        return str(wdl_expr_string.literal)
    parts = wdl_expr_string.parts
    q = cast(str, parts[0])[0]
    string = (
        f"{q}{parts[1]}{q}" if isinstance(parts[1], str) else get_expr(parts[1], ctx)
    )
    if parts[2:-1]:
        string += " + " + " + ".join(
            f"{q}{part}{q}" if isinstance(part, str) else get_expr(part, ctx)
            for part in parts[2:-1]
        )
    return string


def get_expr_name(wdl_expr: WDL.Expr.Ident) -> str:
    """Extract name from WDL expr."""
    return get_input(wdl_expr.name)


def get_expr(wdl_expr: WDL.Expr.Base, ctx: ConversionContext) -> str:
    """Translate WDL Expressions."""
    if isinstance(wdl_expr, WDL.Expr.Apply):
        return get_expr_apply(wdl_expr, ctx)
    elif isinstance(wdl_expr, WDL.Expr.Get):
        return get_expr_get(wdl_expr, ctx)
    elif isinstance(wdl_expr, WDL.Expr.IfThenElse):
        return get_expr_ifthenelse(wdl_expr, ctx)
    elif isinstance(wdl_expr, WDL.Expr.Placeholder):
        return translate_wdl_placeholder(wdl_expr, ctx)
    elif isinstance(wdl_expr, WDL.Expr.String):
        return get_expr_string(wdl_expr, ctx)
    elif isinstance(wdl_expr, WDL.Expr.Boolean) and wdl_expr.literal:
        return str(wdl_expr.literal)  # "true" not "True"
    elif (
        isinstance(
            wdl_expr,
            (
                WDL.Expr.Boolean,
                WDL.Expr.Int,
                WDL.Expr.Float,
            ),
        )
        and wdl_expr.literal
    ):
        return str(wdl_expr.literal.value)
    elif isinstance(wdl_expr, WDL.Expr.Array):
        return "[ " + ", ".join(get_expr(item, ctx) for item in wdl_expr.items) + " ]"
    elif isinstance(wdl_expr, WDL.Expr.Map):
        return (
            "{ "
            + ", ".join(
                f"{get_expr(key, ctx)}: {get_expr(value, ctx)}"
                for key, value in wdl_expr.items
            )
            + " }"
        )
    else:  # pragma: no cover
        raise WDLSourceLine(wdl_expr, ConversionException).makeError(
            f"The expression '{wdl_expr}' is not handled yet."
        )


_BINARY_OPS = {
    "_gt": ">",
    "_lor": "||",
    "_neq": "!==",
    "_lt": "<",
    "_mul": "*",
    "_eqeq": "===",
    "_div": "/",
    "_sub": "-",
}

_SINGLE_ARG_FN = {  # implemented elsewhere, just return the argument
    "read_string",
    "read_float",
    "glob",
    "read_int",
    "read_boolean",
    "read_tsv",
    "read_lines",
}


def get_expr_apply(wdl_apply_expr: WDL.Expr.Apply, ctx: ConversionContext) -> str:
    """Translate WDL Apply Expressions."""
    # N.B: This import here avoids circular dependency error when loading the modules.
    from wdl2cwl import functions

    function_name = wdl_apply_expr.function_name
    arguments = wdl_apply_expr.arguments
    if not arguments:
        raise WDLSourceLine(wdl_apply_expr, ConversionException).makeError(
            f"The '{wdl_apply_expr}' expression has no arguments."
        )
    treat_as_optional = wdl_apply_expr.type.optional

    if function_name in _BINARY_OPS:
        left_operand, right_operand = arguments
        left_operand_expr = get_expr(left_operand, ctx)
        right_operand_expr = get_expr(right_operand, ctx)
        return f"{left_operand_expr} {_BINARY_OPS[function_name]} {right_operand_expr}"
    elif function_name in _SINGLE_ARG_FN:
        only_arg = arguments[0]
        return get_expr(only_arg, ctx)
    elif hasattr(functions, function_name):
        # Call the function if we have it in our wdl2cwl.functions module
        kwargs = {
            "treat_as_optional": treat_as_optional,
            "wdl_apply_expr": wdl_apply_expr,
        }
        return cast(
            str,
            getattr(functions, function_name)(arguments, ctx, **kwargs),
        )
    raise WDLSourceLine(wdl_apply_expr, ConversionException).makeError(
        f"Function name '{function_name}' not yet handled."
    )


def get_expr_get(wdl_get_expr: WDL.Expr.Get, ctx: ConversionContext) -> str:
    """Translate WDL Get Expressions."""
    member = wdl_get_expr.member

    if not member:
        return get_expr_ident(wdl_get_expr.expr, ctx)  # type: ignore[arg-type]
    struct_name = get_expr(wdl_get_expr.expr, ctx)
    member_str = f"{struct_name}.{member}"
    return (
        member_str
        if not isinstance(wdl_get_expr.type, WDL.Type.File)
        else f"{member_str}.path"
    )


def get_expr_ident(wdl_ident_expr: WDL.Expr.Ident, ctx: ConversionContext) -> str:
    """Translate WDL Ident Expressions."""
    id_name = wdl_ident_expr.name
    referee = wdl_ident_expr.referee
    optional = wdl_ident_expr.type.optional
    if referee:
        with WDLSourceLine(referee, ConversionException):
            if isinstance(referee, WDL.Tree.Call):
                return id_name
            if referee.expr and (
                wdl_ident_expr.name in ctx.optional_cwl_null
                or wdl_ident_expr.name not in ctx.non_static_values
            ):
                return get_expr(referee.expr, ctx)
    ident_name = get_input(id_name)
    if optional and isinstance(wdl_ident_expr.type, WDL.Type.File):
        # To prevent null showing on the terminal for inputs of type File
        name_with_file_check = get_expr_name_with_is_file_check(wdl_ident_expr)
        return f'{ident_name} === null ? "" : {name_with_file_check}'
    return (
        ident_name
        if not isinstance(wdl_ident_expr.type, WDL.Type.File)
        else f"{ident_name}.path"
    )


def get_expr_name_with_is_file_check(wdl_expr: WDL.Expr.Ident) -> str:
    """Extract name from WDL expr and check if it's a file path."""
    expr_name = get_input(wdl_expr.name)
    is_file = isinstance(wdl_expr.type, WDL.Type.File)
    return expr_name if not is_file else f"{expr_name}.path"


def get_literal_value(expr: WDL.Expr.Base) -> Optional[Any]:
    """Recursively get a literal value."""
    literal = expr.literal
    if literal:
        if hasattr(expr.parent, "type") and isinstance(expr.parent.type, WDL.Type.File):  # type: ignore[attr-defined]
            return {"class": "File", "path": literal.value}
        value = literal.value
        if isinstance(expr.type, WDL.Type.Map):
            return {key.value: val.value for key, val in value}
        if isinstance(value, list):
            result = []
            for item in value:
                if hasattr(expr.parent, "type") and isinstance(expr.parent.type.item_type, WDL.Type.File):  # type: ignore[attr-defined]
                    result.append({"class": "File", "path": item.value})
                else:
                    result.append(item.value)
            return result
        return value
    return None


def get_step_input_expr(
    wf_expr: Union[WDL.Expr.Get, WDL.Expr.String], ctx: ConversionContext
) -> Tuple[str, Optional[str]]:
    """
    Get name of expression referenced in workflow call inputs.

    Returns a tuple of the source plus any needed "valueFrom" expression.
    """
    with WDLSourceLine(wf_expr, ConversionException):
        if isinstance(wf_expr, WDL.Expr.String):
            return get_expr_string(wf_expr, ctx)[1:-1], None
        elif isinstance(wf_expr, WDL.Expr.Get):
            if isinstance(wf_expr.expr, WDL.Expr.Ident):
                member = None
                id_name = wf_expr.expr.name
                referee = wf_expr.expr.referee
                if referee and isinstance(referee, WDL.Tree.Scatter):
                    scatter_name, value_from = get_step_input_expr(referee.expr, ctx)  # type: ignore[arg-type]
                    ctx.scatter_names.append(scatter_name)
                    return scatter_name, value_from
                return id_name, None
            elif isinstance(wf_expr.expr, WDL.Expr.Get):
                member = str(wf_expr.member)
                ident = cast(WDL.Expr.Ident, wf_expr.expr.expr)
                id_name = ident.name
        elif isinstance(wf_expr, WDL.Expr.Apply):
            expr_str = get_expr(wf_expr, ctx)
            if expr_str.count("inputs") == 1:
                id_name = re.match(r"inputs\.*?[ \.](.*?)[. ]", expr_str).groups()[0]
                value_from = "self" + expr_str.partition(f"inputs.{id_name}")[2]
                return id_name, value_from
        else:
            return get_literal_value(wf_expr), None
        return id_name, f"self.{member}" if member else None


__all__ = [
    "get_expr",
    "get_expr_string",
    "get_step_input_expr",
    "translate_wdl_placeholder",
]
