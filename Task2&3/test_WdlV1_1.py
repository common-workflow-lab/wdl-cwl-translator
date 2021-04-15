f = open("sample.txt", "r")

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

ast = WdlV1_1ParserVisitor().visitDocument(tree)
for a in ast:
    EvaluateTaskVisitor().visit(a)
