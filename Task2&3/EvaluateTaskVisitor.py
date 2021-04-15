
from WdlV1_1Ast import *

class EvaluateTaskVisitor():
    def visit(self, node):
        if type(node) == TaskNode:
            print('Task name: ',node.task_name)
            