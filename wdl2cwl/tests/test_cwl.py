import pytest
import filecmp

import sys
sys.path.append("..")

from .. import WdlV1_1 as wdl

from subprocess import call


class Test:

    def testCollectQualityYieldMetrics(self):
        abc = wdl.main(["wdl_files/Qc.wdl"])
        testStr = ""
        with open('cwl_files/CollectQualityYieldMetrics.cwl','r') as file:
            testStr = file.read()

        assert abc == testStr
