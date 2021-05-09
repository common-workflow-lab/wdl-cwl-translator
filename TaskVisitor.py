from WdlV1_1Ast import *

class TaskVisitor():
    def visit(self, node):
        if type(node) == TaskNode:
            print('Task name: ',node.task_name)
        else:
            print("Node is not Task Node.")

