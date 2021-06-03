f = open("Qc.wdl", "r")

import sys
from antlr4 import *
from WdlV1_1Lexer import WdlV1_1Lexer
from WdlV1_1Parser import WdlV1_1Parser
from WdlV1_1ParserVisitor import WdlV1_1ParserVisitor

text = InputStream(f.read())
lexer = WdlV1_1Lexer(text)
stream = CommonTokenStream(lexer)
parser = WdlV1_1Parser(stream)
tree = parser.document()

#create WdlV1_1ParserVisitor object and return all inputs, outputs, etc
ast = WdlV1_1ParserVisitor()
ast.walk_tree(tree)

#get input type and name
def get_inputs(ast):
    print('TASK INPUTS')
    for a in ast.task_inputs:
        print(a[0], " ", a[1])

#get runtime requirements
def get_runtime(ast):
    print('TASK RUNTIME')
    for a in ast.task_runtime:
        print(a[0]," ",a[1])

#get output type, name, expression
def get_outputs(ast):
    print('TASK OUTPUT')
    for a in ast.task_outputs:
        print(a[0]," ",a[1]," ",a[2])

#not complete. arguments and base command are not separated
def get_command(ast):
    print('TASK COMMAND')
    print(ast.task_command)

def get_arguments(ast):
    print(ast.task_command)

get_inputs(ast)
get_runtime(ast)
get_outputs(ast)
get_command(ast)



