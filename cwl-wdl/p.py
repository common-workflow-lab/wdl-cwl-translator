import sys
from antlr4 import *
from WdlV1_1Lexer import WdlV1_1Lexer
from WdlV1_1Parser import WdlV1_1Parser
from WdlV1_1ParserListener import WdlV1_1ParserListener

class WdlV1_1Listener(WdlV1_1ParserListener):
    def exitWdl_type(self, ctx):
        print("Found")
 
def main(argv):
    lexer = WdlV1_1Lexer(StdinStream)
    stream = CommonTokenStream(lexer)
    parser = WdlV1_1Parser(stream)
    
    printer = WdlV1_1Listener()
    walker = ParseTreeWalker()
    walker.walk(printer, tree)
 
if __name__ == '__main__':
    main(sys.argv)
