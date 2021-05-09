file_document = open("input.txt", "r")

import sys
from antlr4 import *
from WdlV1_1Lexer import WdlV1_1Lexer
from WdlV1_1Parser import WdlV1_1Parser
from WdlV1_1ParserVisitor import WdlV1_1ParserVisitor
from TaskVisitor import TaskVisitor

text = InputStream(file_document.read())
lexer = WdlV1_1Lexer(text)
stream = CommonTokenStream(lexer)
parser = WdlV1_1Parser(stream)
tree = parser.document()
ast = WdlV1_1ParserVisitor().visitDocument(tree)
if ast is None:
    print("Empty AST")
else:
    for a in ast:
        TaskVisitor().visit(a)
