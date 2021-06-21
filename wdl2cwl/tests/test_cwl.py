import os.path


from .. import main as wdl


def get_file(path: str) -> str:
    return os.path.join(os.path.dirname(__file__), path)


def testCollectQualityYieldMetrics() -> None:
    """Test a single WDL to CWL conversion."""
    convertedStr = wdl.main([get_file("wdl_files/CollectQualityYieldMetrics.wdl")])
    testStr = ""
    with open(get_file("cwl_files/CollectQualityYieldMetrics.cwl")) as file:
        testStr = file.read()

    assert convertedStr == testStr

def testSmoove() -> None:
    """Test a single WDL to CWL conversion."""
    convertedStr = wdl.main([get_file("wdl_files/smoove.wdl")])
    testStr = ""
    with open(get_file("cwl_files/smoove.cwl")) as file:
        testStr = file.read()

    assert convertedStr == testStr
