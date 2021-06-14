import pytest
import filecmp

import sys
from subprocess import call
import WdlV1_1 as wdl

import importlib

class Test:

    def testCollectQualityYieldMetrics(self):
        abc = wdl.main()
        testStr = ""
        with open('CollectQualityYieldMetrics.cwl','r') as file:
            testStr = file.read()

        assert abc == testStr
        


#sys.path.append("..")

#from .. import test_WdlV1_1 as wdl
#from wdl-cwl-translator import *

#assert (filecmp.cmp(call('python WdlV1_1.py',  shell=True), "CollectQualityYieldMetrics.cwl", shallow=False))
#assert (filecmp.cmp("../result.cwl", "CollectQualityYieldMetrics.cwl", shallow=False))
#subprocess.call(['python', 'somescript.py'])