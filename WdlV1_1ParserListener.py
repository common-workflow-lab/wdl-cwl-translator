# Generated from WdlV1_1Parser.g4 by ANTLR 4.9.2
from antlr4 import *
if __name__ is not None and "." in __name__:
    from .WdlV1_1Parser import WdlV1_1Parser
else:
    from WdlV1_1Parser import WdlV1_1Parser

# This class defines a complete listener for a parse tree produced by WdlV1_1Parser.
class WdlV1_1ParserListener(ParseTreeListener):

    # Enter a parse tree produced by WdlV1_1Parser#map_type.
    def enterMap_type(self, ctx:WdlV1_1Parser.Map_typeContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#map_type.
    def exitMap_type(self, ctx:WdlV1_1Parser.Map_typeContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#array_type.
    def enterArray_type(self, ctx:WdlV1_1Parser.Array_typeContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#array_type.
    def exitArray_type(self, ctx:WdlV1_1Parser.Array_typeContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#pair_type.
    def enterPair_type(self, ctx:WdlV1_1Parser.Pair_typeContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#pair_type.
    def exitPair_type(self, ctx:WdlV1_1Parser.Pair_typeContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#type_base.
    def enterType_base(self, ctx:WdlV1_1Parser.Type_baseContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#type_base.
    def exitType_base(self, ctx:WdlV1_1Parser.Type_baseContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#wdl_type.
    def enterWdl_type(self, ctx:WdlV1_1Parser.Wdl_typeContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#wdl_type.
    def exitWdl_type(self, ctx:WdlV1_1Parser.Wdl_typeContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#unbound_decls.
    def enterUnbound_decls(self, ctx:WdlV1_1Parser.Unbound_declsContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#unbound_decls.
    def exitUnbound_decls(self, ctx:WdlV1_1Parser.Unbound_declsContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#bound_decls.
    def enterBound_decls(self, ctx:WdlV1_1Parser.Bound_declsContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#bound_decls.
    def exitBound_decls(self, ctx:WdlV1_1Parser.Bound_declsContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#any_decls.
    def enterAny_decls(self, ctx:WdlV1_1Parser.Any_declsContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#any_decls.
    def exitAny_decls(self, ctx:WdlV1_1Parser.Any_declsContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#number.
    def enterNumber(self, ctx:WdlV1_1Parser.NumberContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#number.
    def exitNumber(self, ctx:WdlV1_1Parser.NumberContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#expression_placeholder_option.
    def enterExpression_placeholder_option(self, ctx:WdlV1_1Parser.Expression_placeholder_optionContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#expression_placeholder_option.
    def exitExpression_placeholder_option(self, ctx:WdlV1_1Parser.Expression_placeholder_optionContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#string_part.
    def enterString_part(self, ctx:WdlV1_1Parser.String_partContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#string_part.
    def exitString_part(self, ctx:WdlV1_1Parser.String_partContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#string_expr_part.
    def enterString_expr_part(self, ctx:WdlV1_1Parser.String_expr_partContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#string_expr_part.
    def exitString_expr_part(self, ctx:WdlV1_1Parser.String_expr_partContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#string_expr_with_string_part.
    def enterString_expr_with_string_part(self, ctx:WdlV1_1Parser.String_expr_with_string_partContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#string_expr_with_string_part.
    def exitString_expr_with_string_part(self, ctx:WdlV1_1Parser.String_expr_with_string_partContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#string.
    def enterString(self, ctx:WdlV1_1Parser.StringContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#string.
    def exitString(self, ctx:WdlV1_1Parser.StringContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#primitive_literal.
    def enterPrimitive_literal(self, ctx:WdlV1_1Parser.Primitive_literalContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#primitive_literal.
    def exitPrimitive_literal(self, ctx:WdlV1_1Parser.Primitive_literalContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#expr.
    def enterExpr(self, ctx:WdlV1_1Parser.ExprContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#expr.
    def exitExpr(self, ctx:WdlV1_1Parser.ExprContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#infix0.
    def enterInfix0(self, ctx:WdlV1_1Parser.Infix0Context):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#infix0.
    def exitInfix0(self, ctx:WdlV1_1Parser.Infix0Context):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#infix1.
    def enterInfix1(self, ctx:WdlV1_1Parser.Infix1Context):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#infix1.
    def exitInfix1(self, ctx:WdlV1_1Parser.Infix1Context):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#lor.
    def enterLor(self, ctx:WdlV1_1Parser.LorContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#lor.
    def exitLor(self, ctx:WdlV1_1Parser.LorContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#infix2.
    def enterInfix2(self, ctx:WdlV1_1Parser.Infix2Context):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#infix2.
    def exitInfix2(self, ctx:WdlV1_1Parser.Infix2Context):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#land.
    def enterLand(self, ctx:WdlV1_1Parser.LandContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#land.
    def exitLand(self, ctx:WdlV1_1Parser.LandContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#eqeq.
    def enterEqeq(self, ctx:WdlV1_1Parser.EqeqContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#eqeq.
    def exitEqeq(self, ctx:WdlV1_1Parser.EqeqContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#lt.
    def enterLt(self, ctx:WdlV1_1Parser.LtContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#lt.
    def exitLt(self, ctx:WdlV1_1Parser.LtContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#infix3.
    def enterInfix3(self, ctx:WdlV1_1Parser.Infix3Context):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#infix3.
    def exitInfix3(self, ctx:WdlV1_1Parser.Infix3Context):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#gte.
    def enterGte(self, ctx:WdlV1_1Parser.GteContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#gte.
    def exitGte(self, ctx:WdlV1_1Parser.GteContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#neq.
    def enterNeq(self, ctx:WdlV1_1Parser.NeqContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#neq.
    def exitNeq(self, ctx:WdlV1_1Parser.NeqContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#lte.
    def enterLte(self, ctx:WdlV1_1Parser.LteContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#lte.
    def exitLte(self, ctx:WdlV1_1Parser.LteContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#gt.
    def enterGt(self, ctx:WdlV1_1Parser.GtContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#gt.
    def exitGt(self, ctx:WdlV1_1Parser.GtContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#add.
    def enterAdd(self, ctx:WdlV1_1Parser.AddContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#add.
    def exitAdd(self, ctx:WdlV1_1Parser.AddContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#sub.
    def enterSub(self, ctx:WdlV1_1Parser.SubContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#sub.
    def exitSub(self, ctx:WdlV1_1Parser.SubContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#infix4.
    def enterInfix4(self, ctx:WdlV1_1Parser.Infix4Context):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#infix4.
    def exitInfix4(self, ctx:WdlV1_1Parser.Infix4Context):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#mod.
    def enterMod(self, ctx:WdlV1_1Parser.ModContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#mod.
    def exitMod(self, ctx:WdlV1_1Parser.ModContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#mul.
    def enterMul(self, ctx:WdlV1_1Parser.MulContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#mul.
    def exitMul(self, ctx:WdlV1_1Parser.MulContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#divide.
    def enterDivide(self, ctx:WdlV1_1Parser.DivideContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#divide.
    def exitDivide(self, ctx:WdlV1_1Parser.DivideContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#infix5.
    def enterInfix5(self, ctx:WdlV1_1Parser.Infix5Context):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#infix5.
    def exitInfix5(self, ctx:WdlV1_1Parser.Infix5Context):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#expr_infix5.
    def enterExpr_infix5(self, ctx:WdlV1_1Parser.Expr_infix5Context):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#expr_infix5.
    def exitExpr_infix5(self, ctx:WdlV1_1Parser.Expr_infix5Context):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#member.
    def enterMember(self, ctx:WdlV1_1Parser.MemberContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#member.
    def exitMember(self, ctx:WdlV1_1Parser.MemberContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#pair_literal.
    def enterPair_literal(self, ctx:WdlV1_1Parser.Pair_literalContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#pair_literal.
    def exitPair_literal(self, ctx:WdlV1_1Parser.Pair_literalContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#unarysigned.
    def enterUnarysigned(self, ctx:WdlV1_1Parser.UnarysignedContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#unarysigned.
    def exitUnarysigned(self, ctx:WdlV1_1Parser.UnarysignedContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#apply.
    def enterApply(self, ctx:WdlV1_1Parser.ApplyContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#apply.
    def exitApply(self, ctx:WdlV1_1Parser.ApplyContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#expression_group.
    def enterExpression_group(self, ctx:WdlV1_1Parser.Expression_groupContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#expression_group.
    def exitExpression_group(self, ctx:WdlV1_1Parser.Expression_groupContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#primitives.
    def enterPrimitives(self, ctx:WdlV1_1Parser.PrimitivesContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#primitives.
    def exitPrimitives(self, ctx:WdlV1_1Parser.PrimitivesContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#left_name.
    def enterLeft_name(self, ctx:WdlV1_1Parser.Left_nameContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#left_name.
    def exitLeft_name(self, ctx:WdlV1_1Parser.Left_nameContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#at.
    def enterAt(self, ctx:WdlV1_1Parser.AtContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#at.
    def exitAt(self, ctx:WdlV1_1Parser.AtContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#negate.
    def enterNegate(self, ctx:WdlV1_1Parser.NegateContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#negate.
    def exitNegate(self, ctx:WdlV1_1Parser.NegateContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#map_literal.
    def enterMap_literal(self, ctx:WdlV1_1Parser.Map_literalContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#map_literal.
    def exitMap_literal(self, ctx:WdlV1_1Parser.Map_literalContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#ifthenelse.
    def enterIfthenelse(self, ctx:WdlV1_1Parser.IfthenelseContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#ifthenelse.
    def exitIfthenelse(self, ctx:WdlV1_1Parser.IfthenelseContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#get_name.
    def enterGet_name(self, ctx:WdlV1_1Parser.Get_nameContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#get_name.
    def exitGet_name(self, ctx:WdlV1_1Parser.Get_nameContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#object_literal.
    def enterObject_literal(self, ctx:WdlV1_1Parser.Object_literalContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#object_literal.
    def exitObject_literal(self, ctx:WdlV1_1Parser.Object_literalContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#array_literal.
    def enterArray_literal(self, ctx:WdlV1_1Parser.Array_literalContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#array_literal.
    def exitArray_literal(self, ctx:WdlV1_1Parser.Array_literalContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#struct_literal.
    def enterStruct_literal(self, ctx:WdlV1_1Parser.Struct_literalContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#struct_literal.
    def exitStruct_literal(self, ctx:WdlV1_1Parser.Struct_literalContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#version.
    def enterVersion(self, ctx:WdlV1_1Parser.VersionContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#version.
    def exitVersion(self, ctx:WdlV1_1Parser.VersionContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#import_alias.
    def enterImport_alias(self, ctx:WdlV1_1Parser.Import_aliasContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#import_alias.
    def exitImport_alias(self, ctx:WdlV1_1Parser.Import_aliasContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#import_as.
    def enterImport_as(self, ctx:WdlV1_1Parser.Import_asContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#import_as.
    def exitImport_as(self, ctx:WdlV1_1Parser.Import_asContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#import_doc.
    def enterImport_doc(self, ctx:WdlV1_1Parser.Import_docContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#import_doc.
    def exitImport_doc(self, ctx:WdlV1_1Parser.Import_docContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#struct.
    def enterStruct(self, ctx:WdlV1_1Parser.StructContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#struct.
    def exitStruct(self, ctx:WdlV1_1Parser.StructContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#meta_value.
    def enterMeta_value(self, ctx:WdlV1_1Parser.Meta_valueContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#meta_value.
    def exitMeta_value(self, ctx:WdlV1_1Parser.Meta_valueContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#meta_string_part.
    def enterMeta_string_part(self, ctx:WdlV1_1Parser.Meta_string_partContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#meta_string_part.
    def exitMeta_string_part(self, ctx:WdlV1_1Parser.Meta_string_partContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#meta_string.
    def enterMeta_string(self, ctx:WdlV1_1Parser.Meta_stringContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#meta_string.
    def exitMeta_string(self, ctx:WdlV1_1Parser.Meta_stringContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#meta_array.
    def enterMeta_array(self, ctx:WdlV1_1Parser.Meta_arrayContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#meta_array.
    def exitMeta_array(self, ctx:WdlV1_1Parser.Meta_arrayContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#meta_object.
    def enterMeta_object(self, ctx:WdlV1_1Parser.Meta_objectContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#meta_object.
    def exitMeta_object(self, ctx:WdlV1_1Parser.Meta_objectContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#meta_object_kv.
    def enterMeta_object_kv(self, ctx:WdlV1_1Parser.Meta_object_kvContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#meta_object_kv.
    def exitMeta_object_kv(self, ctx:WdlV1_1Parser.Meta_object_kvContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#meta_kv.
    def enterMeta_kv(self, ctx:WdlV1_1Parser.Meta_kvContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#meta_kv.
    def exitMeta_kv(self, ctx:WdlV1_1Parser.Meta_kvContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#parameter_meta.
    def enterParameter_meta(self, ctx:WdlV1_1Parser.Parameter_metaContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#parameter_meta.
    def exitParameter_meta(self, ctx:WdlV1_1Parser.Parameter_metaContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#meta.
    def enterMeta(self, ctx:WdlV1_1Parser.MetaContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#meta.
    def exitMeta(self, ctx:WdlV1_1Parser.MetaContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#task_runtime_kv.
    def enterTask_runtime_kv(self, ctx:WdlV1_1Parser.Task_runtime_kvContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#task_runtime_kv.
    def exitTask_runtime_kv(self, ctx:WdlV1_1Parser.Task_runtime_kvContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#task_runtime.
    def enterTask_runtime(self, ctx:WdlV1_1Parser.Task_runtimeContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#task_runtime.
    def exitTask_runtime(self, ctx:WdlV1_1Parser.Task_runtimeContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#task_input.
    def enterTask_input(self, ctx:WdlV1_1Parser.Task_inputContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#task_input.
    def exitTask_input(self, ctx:WdlV1_1Parser.Task_inputContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#task_output.
    def enterTask_output(self, ctx:WdlV1_1Parser.Task_outputContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#task_output.
    def exitTask_output(self, ctx:WdlV1_1Parser.Task_outputContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#task_command_string_part.
    def enterTask_command_string_part(self, ctx:WdlV1_1Parser.Task_command_string_partContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#task_command_string_part.
    def exitTask_command_string_part(self, ctx:WdlV1_1Parser.Task_command_string_partContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#task_command_expr_part.
    def enterTask_command_expr_part(self, ctx:WdlV1_1Parser.Task_command_expr_partContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#task_command_expr_part.
    def exitTask_command_expr_part(self, ctx:WdlV1_1Parser.Task_command_expr_partContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#task_command_expr_with_string.
    def enterTask_command_expr_with_string(self, ctx:WdlV1_1Parser.Task_command_expr_with_stringContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#task_command_expr_with_string.
    def exitTask_command_expr_with_string(self, ctx:WdlV1_1Parser.Task_command_expr_with_stringContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#task_command.
    def enterTask_command(self, ctx:WdlV1_1Parser.Task_commandContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#task_command.
    def exitTask_command(self, ctx:WdlV1_1Parser.Task_commandContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#task_element.
    def enterTask_element(self, ctx:WdlV1_1Parser.Task_elementContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#task_element.
    def exitTask_element(self, ctx:WdlV1_1Parser.Task_elementContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#task.
    def enterTask(self, ctx:WdlV1_1Parser.TaskContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#task.
    def exitTask(self, ctx:WdlV1_1Parser.TaskContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#inner_workflow_element.
    def enterInner_workflow_element(self, ctx:WdlV1_1Parser.Inner_workflow_elementContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#inner_workflow_element.
    def exitInner_workflow_element(self, ctx:WdlV1_1Parser.Inner_workflow_elementContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#call_alias.
    def enterCall_alias(self, ctx:WdlV1_1Parser.Call_aliasContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#call_alias.
    def exitCall_alias(self, ctx:WdlV1_1Parser.Call_aliasContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#call_input.
    def enterCall_input(self, ctx:WdlV1_1Parser.Call_inputContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#call_input.
    def exitCall_input(self, ctx:WdlV1_1Parser.Call_inputContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#call_inputs.
    def enterCall_inputs(self, ctx:WdlV1_1Parser.Call_inputsContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#call_inputs.
    def exitCall_inputs(self, ctx:WdlV1_1Parser.Call_inputsContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#call_body.
    def enterCall_body(self, ctx:WdlV1_1Parser.Call_bodyContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#call_body.
    def exitCall_body(self, ctx:WdlV1_1Parser.Call_bodyContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#call_after.
    def enterCall_after(self, ctx:WdlV1_1Parser.Call_afterContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#call_after.
    def exitCall_after(self, ctx:WdlV1_1Parser.Call_afterContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#call_name.
    def enterCall_name(self, ctx:WdlV1_1Parser.Call_nameContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#call_name.
    def exitCall_name(self, ctx:WdlV1_1Parser.Call_nameContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#call.
    def enterCall(self, ctx:WdlV1_1Parser.CallContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#call.
    def exitCall(self, ctx:WdlV1_1Parser.CallContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#scatter.
    def enterScatter(self, ctx:WdlV1_1Parser.ScatterContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#scatter.
    def exitScatter(self, ctx:WdlV1_1Parser.ScatterContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#conditional.
    def enterConditional(self, ctx:WdlV1_1Parser.ConditionalContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#conditional.
    def exitConditional(self, ctx:WdlV1_1Parser.ConditionalContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#workflow_input.
    def enterWorkflow_input(self, ctx:WdlV1_1Parser.Workflow_inputContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#workflow_input.
    def exitWorkflow_input(self, ctx:WdlV1_1Parser.Workflow_inputContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#workflow_output.
    def enterWorkflow_output(self, ctx:WdlV1_1Parser.Workflow_outputContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#workflow_output.
    def exitWorkflow_output(self, ctx:WdlV1_1Parser.Workflow_outputContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#input.
    def enterInput(self, ctx:WdlV1_1Parser.InputContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#input.
    def exitInput(self, ctx:WdlV1_1Parser.InputContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#output.
    def enterOutput(self, ctx:WdlV1_1Parser.OutputContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#output.
    def exitOutput(self, ctx:WdlV1_1Parser.OutputContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#inner_element.
    def enterInner_element(self, ctx:WdlV1_1Parser.Inner_elementContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#inner_element.
    def exitInner_element(self, ctx:WdlV1_1Parser.Inner_elementContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#parameter_meta_element.
    def enterParameter_meta_element(self, ctx:WdlV1_1Parser.Parameter_meta_elementContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#parameter_meta_element.
    def exitParameter_meta_element(self, ctx:WdlV1_1Parser.Parameter_meta_elementContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#meta_element.
    def enterMeta_element(self, ctx:WdlV1_1Parser.Meta_elementContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#meta_element.
    def exitMeta_element(self, ctx:WdlV1_1Parser.Meta_elementContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#workflow.
    def enterWorkflow(self, ctx:WdlV1_1Parser.WorkflowContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#workflow.
    def exitWorkflow(self, ctx:WdlV1_1Parser.WorkflowContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#document_element.
    def enterDocument_element(self, ctx:WdlV1_1Parser.Document_elementContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#document_element.
    def exitDocument_element(self, ctx:WdlV1_1Parser.Document_elementContext):
        pass


    # Enter a parse tree produced by WdlV1_1Parser#document.
    def enterDocument(self, ctx:WdlV1_1Parser.DocumentContext):
        pass

    # Exit a parse tree produced by WdlV1_1Parser#document.
    def exitDocument(self, ctx:WdlV1_1Parser.DocumentContext):
        pass



del WdlV1_1Parser