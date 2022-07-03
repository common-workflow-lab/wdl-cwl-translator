"""WDL Stdlib Functions."""

import re
from typing import List, Any, cast

import WDL
from WDL.Expr import Base
from wdl2cwl.errors import WDLSourceLine, ConversionException
from wdl2cwl.expr import (
    get_expr,
    get_expr_get,
    get_expr_name,
    get_expr_name_with_is_file_check,
    get_literal_value,
)
from wdl2cwl.util import ConversionContext, get_input, get_mem_in_bytes


def _add(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib add function."""
    add_left_operand = arguments[0]
    add_right_operand = get_expr(arguments[1], ctx)
    add_left_operand_value = get_expr(add_left_operand, ctx)
    treat_as_optional = (
        kwargs["treat_as_optional"] if "treat_as_optional" in kwargs else False
    )
    referer = ""
    if getattr(add_left_operand, "function_name", None) == "basename":
        referer = kwargs["wdl_apply_expr"].parent.name
        treat_as_optional = True if referer in ctx.non_static_values else False
    return (
        f"{add_left_operand_value} + {add_right_operand}"
        if not treat_as_optional
        else f"{get_input(referer)} === null ? {add_left_operand_value} + {add_right_operand} : {get_input(referer)}"
    )


def basename(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib basename function."""
    if len(arguments) == 1:
        only_operand = arguments[0]
        is_file = isinstance(only_operand.type, WDL.Type.File)
        if isinstance(only_operand, WDL.Expr.Get) and isinstance(
            only_operand.expr, WDL.Expr.Ident
        ):
            only_operand_name = get_expr_name(only_operand.expr)
        else:
            only_operand_name = get_expr(only_operand, ctx)
        return (
            f"{only_operand_name}.basename"
            if is_file
            else f"{only_operand_name}.split('/').reverse()[0]"
        )
    else:
        basename_target, suffix = arguments
        is_file = isinstance(basename_target.type, WDL.Type.File)
        if isinstance(basename_target, WDL.Expr.Get):
            basename_target_name = get_expr_name(basename_target.expr)  # type: ignore[arg-type]
        elif isinstance(basename_target, WDL.Expr.Apply):
            basename_target_name = get_expr(basename_target, ctx)
        suffix_str = str(get_literal_value(suffix))
        regex_str = re.escape(suffix_str)
        return (
            f"{basename_target_name}.basename.replace(/{regex_str}$/, '') "
            if is_file
            else f"{basename_target_name}.split('/').reverse()[0].replace(/{regex_str}$/, '')"
        )


def defined(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib defined function."""
    only_operand = arguments[0]
    assert isinstance(only_operand, WDL.Expr.Get) and isinstance(
        only_operand.expr, WDL.Expr.Ident
    )
    return f"{get_expr_name(only_operand.expr)} !== null"


def _interpolation_add(
    arguments: List[Base], ctx: ConversionContext, **kwargs: Any
) -> str:
    """WDL Stdlib _interpolation_add function."""
    arg_value, arg_name = arguments
    if isinstance(arg_name, WDL.Expr.String) and isinstance(
        arg_value, (WDL.Expr.Apply, WDL.Expr.String)
    ):
        return f"{get_expr(arg_value, ctx)} + {get_expr(arg_name, ctx)}"
    treat_as_optional = (
        kwargs["treat_as_optional"] if "treat_as_optional" in kwargs else False
    )
    if isinstance(arg_name, (WDL.Expr.Placeholder, WDL.Expr.Get)):
        just_arg_name = get_expr_name(arg_name.expr)  # type: ignore[arg-type]
        arg_name_with_file_check = get_expr_name_with_is_file_check(
            arg_name.expr  # type: ignore[arg-type]
        )
    elif isinstance(arg_value, (WDL.Expr.Placeholder, WDL.Expr.Get)):
        just_arg_name = get_expr_name(arg_value.expr)  # type: ignore[arg-type]
        arg_name_with_file_check = get_expr_name_with_is_file_check(
            arg_value.expr  # type: ignore[arg-type]
        )
        arg_value = arg_name
    with WDLSourceLine(arg_value, ConversionException):
        arg_value_str = get_expr(arg_value, ctx)
        return (
            f'{just_arg_name} === null ? "" : {arg_value_str} + {arg_name_with_file_check}'
            if treat_as_optional
            else f"{arg_value_str} + {arg_name_with_file_check}"
        )


def sub(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib sub function."""
    wdl_apply, arg_string, arg_sub = arguments
    sub_expr = get_expr(wdl_apply, ctx)
    arg_string_expr = get_expr(arg_string, ctx)
    arg_sub_expr = get_expr(arg_sub, ctx)
    return f"{sub_expr}.replace({arg_string_expr}, {arg_sub_expr}) "


def _at(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib _at function."""
    iterable_object, index = arguments
    iterable_object_expr = get_expr(iterable_object, ctx)
    index_expr = get_expr(index, ctx)
    return f"{iterable_object_expr}[{index_expr}]"


def length(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib length function."""
    only_arg_expr = get_expr_get(arguments[0], ctx)  # type: ignore[arg-type]
    return f"{only_arg_expr}.length"


def round(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib round function."""
    only_arg_expr = get_expr(arguments[0], ctx)
    return f"Math.round({only_arg_expr})"


def select_first(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib select_first function."""
    array_obj = cast(WDL.Expr.Array, arguments[0])
    array_items = [str(get_expr(item, ctx)) for item in array_obj.items]
    items_str = ", ".join(array_items)
    return f"[{items_str}].find(function(element) {{ return element !== null }}) "


def select_all(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib select_all function."""
    array_obj = cast(WDL.Expr.Array, arguments[0])
    array_items = [str(get_expr(item, ctx)) for item in array_obj.items]
    items_str = ", ".join(array_items)
    return f"[{items_str}].filter(function(element) {{ return element !== null }}) "


def ceil(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib ceil function."""
    only_arg = get_expr(arguments[0], ctx)
    return f"Math.ceil({only_arg}) "


def size(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib size function."""
    if len(arguments) == 1:
        left_operand = arguments[0]
        unit_value = "1"
    else:
        left_operand, right_operand = arguments
        right_value = str(get_literal_value(right_operand))
        unit_base, unit_exponent = get_mem_in_bytes(right_value)
        unit_value = f"{unit_base}^{unit_exponent}"
    if isinstance(left_operand, WDL.Expr.Array):
        array_items = [get_expr(item, ctx) for item in left_operand.items]
        left = ", ".join(array_items)
        left_str = f"[{left}]"
    else:
        left_str = get_expr(left_operand, ctx)
    return (
        "(function(size_of=0)"
        + "{"
        + f"{left_str}.forEach(function(element)"
        + "{ if (element) {"
        + "size_of += element.size"
        + "}})}"
        + f") / {unit_value}"
    )


def flatten(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib flatten function."""
    flatten_array = arguments[0]
    with WDLSourceLine(flatten_array, ConversionException):
        items_str = get_expr(flatten_array, ctx)
    result = (
        "(function () {var new_array = []; "
        + items_str
        + ".forEach(function(value, index, obj) "
        "{value.forEach(function(sub_value, sub_index, sub_obj) "
        "{new_array.push(sub_value);});}); return new_array;})()"
    )
    return result


def sep(arguments: List[Base], ctx: ConversionContext, **kwargs: Any) -> str:
    """WDL Stdlib sep function."""
    sep, array = arguments
    if isinstance(array, WDL.Expr.Get) and isinstance(array.type, WDL.Type.Array):
        item_type = array.type.item_type
    else:
        raise WDLSourceLine(array, ConversionException).makeError(
            f"Unhandled sep array type: {type(array)}: {array}."
        )
    sep_str = get_literal_value(sep) or ""
    if isinstance(item_type, WDL.Type.File):
        return (
            f"{get_expr(array, ctx)}.map("
            + 'function(el) {return el.path}).join("'
            + sep_str
            + '")'
        )
    else:
        return f'{get_expr(array, ctx)}.join("{sep_str}")'


__all__ = [
    "_add",
    "basename",
    "defined",
    "_interpolation_add",
    "sub",
    "_at",
    "length",
    "round",
    "select_first",
    "select_all",
    "ceil",
    "size",
    "flatten",
    "sep",
]
