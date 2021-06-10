import pytest

import filecmp

class Test:

    def testCollectQualityYieldMetrics(self):
        assert (filecmp.cmp("../result.cwl", "CollectQualityYieldMetrics.cwl", shallow=False))