# Generated from WdlV1_1Parser.g4 by ANTLR 4.7.2
from antlr4 import *
if __name__ is not None and "." in __name__:
    from .WdlV1_1Parser import WdlV1_1Parser
else:
    from WdlV1_1Parser import WdlV1_1Parser

# This class defines a complete generic visitor for a parse tree produced by WdlV1_1Parser.

class WdlV1_1ParserVisitor(ParseTreeVisitor):

    def __init__(self):
        self.task_inputs = []
        self.task_outputs = []
        self.task_command = None
        self.task_runtime = []
        #checks are used to check the parent node
        self.task_input_check = None
        self.workflow_input_check = None
        self.task_output_check = None

    def walk_tree(self, tree):
        self.visitDocument(tree)
        return self

    # Visit a parse tree produced by WdlV1_1Parser#map_type.
    def visitMap_type(self, ctx:WdlV1_1Parser.Map_typeContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#array_type.
    def visitArray_type(self, ctx:WdlV1_1Parser.Array_typeContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#pair_type.
    def visitPair_type(self, ctx:WdlV1_1Parser.Pair_typeContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#type_base.
    def visitType_base(self, ctx:WdlV1_1Parser.Type_baseContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#wdl_type.
    def visitWdl_type(self, ctx:WdlV1_1Parser.Wdl_typeContext):
        return ctx.getText()
        #return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#unbound_decls.
    def visitUnbound_decls(self, ctx:WdlV1_1Parser.Unbound_declsContext):
        input_type = self.visitWdl_type(ctx.wdl_type())
        return [input_type, ctx.Identifier()]

    # Visit a parse tree produced by WdlV1_1Parser#bound_decls.
    def visitBound_decls(self, ctx:WdlV1_1Parser.Bound_declsContext):
        decl_type = self.visitWdl_type(ctx.wdl_type())
        expression = self.visitExpr(ctx.expr())
        if self.task_output_check:
            self.task_outputs.append([decl_type, ctx.Identifier(), expression])
        return self.visitChildren(ctx)
        
    # Visit a parse tree produced by WdlV1_1Parser#any_decls.
    def visitAny_decls(self, ctx:WdlV1_1Parser.Any_declsContext):
        unbound_decls = self.visitUnbound_decls(ctx.unbound_decls())
        if self.task_input_check:
            self.task_inputs.append(unbound_decls)

        return unbound_decls

    # Visit a parse tree produced by WdlV1_1Parser#number.
    def visitNumber(self, ctx:WdlV1_1Parser.NumberContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#expression_placeholder_option.
    def visitExpression_placeholder_option(self, ctx:WdlV1_1Parser.Expression_placeholder_optionContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#string_part.
    def visitString_part(self, ctx:WdlV1_1Parser.String_partContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#string_expr_part.
    def visitString_expr_part(self, ctx:WdlV1_1Parser.String_expr_partContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#string_expr_with_string_part.
    def visitString_expr_with_string_part(self, ctx:WdlV1_1Parser.String_expr_with_string_partContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#string.
    def visitString(self, ctx:WdlV1_1Parser.StringContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#primitive_literal.
    def visitPrimitive_literal(self, ctx:WdlV1_1Parser.Primitive_literalContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#expr.
    def visitExpr(self, ctx:WdlV1_1Parser.ExprContext):
        return ctx.getText()


    # Visit a parse tree produced by WdlV1_1Parser#infix0.
    def visitInfix0(self, ctx:WdlV1_1Parser.Infix0Context):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#infix1.
    def visitInfix1(self, ctx:WdlV1_1Parser.Infix1Context):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#lor.
    def visitLor(self, ctx:WdlV1_1Parser.LorContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#infix2.
    def visitInfix2(self, ctx:WdlV1_1Parser.Infix2Context):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#land.
    def visitLand(self, ctx:WdlV1_1Parser.LandContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#eqeq.
    def visitEqeq(self, ctx:WdlV1_1Parser.EqeqContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#lt.
    def visitLt(self, ctx:WdlV1_1Parser.LtContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#infix3.
    def visitInfix3(self, ctx:WdlV1_1Parser.Infix3Context):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#gte.
    def visitGte(self, ctx:WdlV1_1Parser.GteContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#neq.
    def visitNeq(self, ctx:WdlV1_1Parser.NeqContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#lte.
    def visitLte(self, ctx:WdlV1_1Parser.LteContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#gt.
    def visitGt(self, ctx:WdlV1_1Parser.GtContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#add.
    def visitAdd(self, ctx:WdlV1_1Parser.AddContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#sub.
    def visitSub(self, ctx:WdlV1_1Parser.SubContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#infix4.
    def visitInfix4(self, ctx:WdlV1_1Parser.Infix4Context):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#mod.
    def visitMod(self, ctx:WdlV1_1Parser.ModContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#mul.
    def visitMul(self, ctx:WdlV1_1Parser.MulContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#divide.
    def visitDivide(self, ctx:WdlV1_1Parser.DivideContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#infix5.
    def visitInfix5(self, ctx:WdlV1_1Parser.Infix5Context):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#expr_infix5.
    def visitExpr_infix5(self, ctx:WdlV1_1Parser.Expr_infix5Context):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#member.
    def visitMember(self, ctx:WdlV1_1Parser.MemberContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#pair_literal.
    def visitPair_literal(self, ctx:WdlV1_1Parser.Pair_literalContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#unarysigned.
    def visitUnarysigned(self, ctx:WdlV1_1Parser.UnarysignedContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#apply.
    def visitApply(self, ctx:WdlV1_1Parser.ApplyContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#expression_group.
    def visitExpression_group(self, ctx:WdlV1_1Parser.Expression_groupContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#primitives.
    def visitPrimitives(self, ctx:WdlV1_1Parser.PrimitivesContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#left_name.
    def visitLeft_name(self, ctx:WdlV1_1Parser.Left_nameContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#at.
    def visitAt(self, ctx:WdlV1_1Parser.AtContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#negate.
    def visitNegate(self, ctx:WdlV1_1Parser.NegateContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#map_literal.
    def visitMap_literal(self, ctx:WdlV1_1Parser.Map_literalContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#ifthenelse.
    def visitIfthenelse(self, ctx:WdlV1_1Parser.IfthenelseContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#get_name.
    def visitGet_name(self, ctx:WdlV1_1Parser.Get_nameContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#object_literal.
    def visitObject_literal(self, ctx:WdlV1_1Parser.Object_literalContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#array_literal.
    def visitArray_literal(self, ctx:WdlV1_1Parser.Array_literalContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#struct_literal.
    def visitStruct_literal(self, ctx:WdlV1_1Parser.Struct_literalContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#version.
    def visitVersion(self, ctx:WdlV1_1Parser.VersionContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#import_alias.
    def visitImport_alias(self, ctx:WdlV1_1Parser.Import_aliasContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#import_as.
    def visitImport_as(self, ctx:WdlV1_1Parser.Import_asContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#import_doc.
    def visitImport_doc(self, ctx:WdlV1_1Parser.Import_docContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#struct.
    def visitStruct(self, ctx:WdlV1_1Parser.StructContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#meta_value.
    def visitMeta_value(self, ctx:WdlV1_1Parser.Meta_valueContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#meta_string_part.
    def visitMeta_string_part(self, ctx:WdlV1_1Parser.Meta_string_partContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#meta_string.
    def visitMeta_string(self, ctx:WdlV1_1Parser.Meta_stringContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#meta_array.
    def visitMeta_array(self, ctx:WdlV1_1Parser.Meta_arrayContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#meta_object.
    def visitMeta_object(self, ctx:WdlV1_1Parser.Meta_objectContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#meta_object_kv.
    def visitMeta_object_kv(self, ctx:WdlV1_1Parser.Meta_object_kvContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#meta_kv.
    def visitMeta_kv(self, ctx:WdlV1_1Parser.Meta_kvContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#parameter_meta.
    def visitParameter_meta(self, ctx:WdlV1_1Parser.Parameter_metaContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#meta.
    def visitMeta(self, ctx:WdlV1_1Parser.MetaContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#task_runtime_kv.
    def visitTask_runtime_kv(self, ctx:WdlV1_1Parser.Task_runtime_kvContext):
        expression = self.visitExpr(ctx.expr())      
        self.task_runtime.append([ctx.Identifier(), expression])
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#task_runtime.
    def visitTask_runtime(self, ctx:WdlV1_1Parser.Task_runtimeContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#task_input.
    def visitTask_input(self, ctx:WdlV1_1Parser.Task_inputContext):
        self.task_input_check = 1
        
        return self.visitChildren(ctx)
        
    # Visit a parse tree produced by WdlV1_1Parser#task_output.
    def visitTask_output(self, ctx:WdlV1_1Parser.Task_outputContext):
        self.task_output_check = 1
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#task_command_string_part.
    def visitTask_command_string_part(self, ctx:WdlV1_1Parser.Task_command_string_partContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#task_command_expr_part.
    def visitTask_command_expr_part(self, ctx:WdlV1_1Parser.Task_command_expr_partContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#task_command_expr_with_string.
    def visitTask_command_expr_with_string(self, ctx:WdlV1_1Parser.Task_command_expr_with_stringContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#task_command.
    def visitTask_command(self, ctx:WdlV1_1Parser.Task_commandContext):
        self.task_command=(ctx.getText())
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#task_element.
    def visitTask_element(self, ctx:WdlV1_1Parser.Task_elementContext):   
        self.task_input_check = None
        self.task_output_check = None
        return self.visitChildren(ctx)

    # Visit a parse tree produced by WdlV1_1Parser#task.
    def visitTask(self, ctx:WdlV1_1Parser.TaskContext):
        return self.visitChildren(ctx)

    # Visit a parse tree produced by WdlV1_1Parser#inner_workflow_element.
    def visitInner_workflow_element(self, ctx:WdlV1_1Parser.Inner_workflow_elementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#call_alias.
    def visitCall_alias(self, ctx:WdlV1_1Parser.Call_aliasContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#call_input.
    def visitCall_input(self, ctx:WdlV1_1Parser.Call_inputContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#call_inputs.
    def visitCall_inputs(self, ctx:WdlV1_1Parser.Call_inputsContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#call_body.
    def visitCall_body(self, ctx:WdlV1_1Parser.Call_bodyContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#call_after.
    def visitCall_after(self, ctx:WdlV1_1Parser.Call_afterContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#call_name.
    def visitCall_name(self, ctx:WdlV1_1Parser.Call_nameContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#call.
    def visitCall(self, ctx:WdlV1_1Parser.CallContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#scatter.
    def visitScatter(self, ctx:WdlV1_1Parser.ScatterContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#conditional.
    def visitConditional(self, ctx:WdlV1_1Parser.ConditionalContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#workflow_input.
    def visitWorkflow_input(self, ctx:WdlV1_1Parser.Workflow_inputContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#workflow_output.
    def visitWorkflow_output(self, ctx:WdlV1_1Parser.Workflow_outputContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#input.
    def visitInput(self, ctx:WdlV1_1Parser.InputContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#output.
    def visitOutput(self, ctx:WdlV1_1Parser.OutputContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#inner_element.
    def visitInner_element(self, ctx:WdlV1_1Parser.Inner_elementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#parameter_meta_element.
    def visitParameter_meta_element(self, ctx:WdlV1_1Parser.Parameter_meta_elementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#meta_element.
    def visitMeta_element(self, ctx:WdlV1_1Parser.Meta_elementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by WdlV1_1Parser#workflow.
    def visitWorkflow(self, ctx:WdlV1_1Parser.WorkflowContext):
        return self.visitChildren(ctx)

    # Visit a parse tree produced by WdlV1_1Parser#document_element.
    def visitDocument_element(self, ctx:WdlV1_1Parser.Document_elementContext):
        return self.visitChildren(ctx)

    # Visit a parse tree produced by WdlV1_1Parser#document.
    def visitDocument(self, ctx:WdlV1_1Parser.DocumentContext):
        return self.visitChildren(ctx)
        

del WdlV1_1Parser
