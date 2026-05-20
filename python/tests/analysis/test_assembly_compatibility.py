import pandas as pd

from bsx2.analysis.assembly_compatibility import check_coordinate_bounds


def test_valid_coordinates_pass():
    report = check_coordinate_bounds(pd.DataFrame({"seqname": ["A"], "start": [0], "end": [5]}), pd.DataFrame({"seqname": ["A"], "length": [10]}))
    assert report["compatibility_status"].iloc[0] == "ok"


def test_out_of_bounds_flagged():
    report = check_coordinate_bounds(pd.DataFrame({"seqname": ["A"], "start": [0], "end": [50]}), pd.DataFrame({"seqname": ["A"], "length": [10]}))
    assert report["out_of_bounds"].iloc[0]
