f = open("Qc.wdl", "r")

import sys
from antlr4 import *
from WdlV1_1Lexer import WdlV1_1Lexer
from WdlV1_1Parser import WdlV1_1Parser
from WdlV1_1ParserVisitor import WdlV1_1ParserVisitor
from EvaluateTaskVisitor import EvaluateTaskVisitor

text = InputStream(f.read())
lexer = WdlV1_1Lexer(text)
stream = CommonTokenStream(lexer)
parser = WdlV1_1Parser(stream)
tree = parser.document()
tree2 = parser.task_input()
tree3 = parser.unbound_decls()

def visit_tasks():
    ast = WdlV1_1ParserVisitor().visitDocument(tree)
    for a in ast:
        EvaluateTaskVisitor().visit(a)

def get_command():
    
    ast = WdlV1_1ParserVisitor().visitDocument(tree)

get_command()
#visit_tasks()
#get_input()

