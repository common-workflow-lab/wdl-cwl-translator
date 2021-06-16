import os.path


from .. import main as wdl


def get_file(path: str) -> str:
    return os.path.join(os.path.dirname(__file__), path)


def testCollectQualityYieldMetrics() -> None:
    """Test a single WDL to CWL conversion."""
    abc = wdl.main([get_file("wdl_files/Qc.wdl")])
    testStr = ""
    with open(get_file("cwl_files/CollectQualityYieldMetrics.cwl")) as file:
        testStr = file.read()

    assert abc == testStr
