import pytest
import filecmp

import os.path
import sys

from .. import WdlV1_1 as wdl

from subprocess import call

def get_file(path: str) -> str:
    return os.path.join(os.path.dirname(__file__), path)

class Test:

    def testCollectQualityYieldMetrics(self):
        abc = wdl.main([get_file("wdl_files/Qc.wdl")])
        testStr = ""
        with open(get_file('cwl_files/CollectQualityYieldMetrics.cwl'),'r') as file:
            testStr = file.read()

        assert abc == testStr
